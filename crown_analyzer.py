#!/usr/bin/env python3
"""
Crown Adjustment Analysis Tool v1.3
===================================
- Adjustable brightness control
- Vertical scalar bar (GeoMagic style)
- Even lighting across all surfaces
"""

import sys
import os
import numpy as np
from pathlib import Path
from datetime import datetime

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QGroupBox, QGridLayout,
    QTabWidget, QTextEdit, QProgressBar, QMessageBox, QFrame,
    QSplitter, QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox,
    QSlider
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont

import trimesh
from scipy.spatial import cKDTree
import pyvista as pv
from pyvistaqt import QtInteractor
from typing import Optional, Tuple, Dict


class MeshAnalyzer:
    def __init__(self):
        self.baseline_mesh: Optional[trimesh.Trimesh] = None
        self.adjusted_mesh: Optional[trimesh.Trimesh] = None
        self.distances: Optional[np.ndarray] = None
        self.closest_points: Optional[np.ndarray] = None
        
    def load_mesh(self, filepath: str) -> Tuple[trimesh.Trimesh, dict]:
        mesh = trimesh.load(filepath, force='mesh')
        if not mesh.is_watertight:
            mesh.fill_holes()
            mesh.fix_normals()
        info = {
            'filepath': filepath,
            'filename': os.path.basename(filepath),
            'vertices': len(mesh.vertices),
            'faces': len(mesh.faces),
            'is_watertight': mesh.is_watertight,
            'volume_mm3': abs(mesh.volume) if mesh.is_watertight else None,
            'surface_area_mm2': mesh.area,
            'bounds': mesh.bounds,
            'centroid': mesh.centroid
        }
        return mesh, info
    
    def calculate_volume(self, mesh: trimesh.Trimesh) -> float:
        if mesh.is_watertight:
            return abs(mesh.volume)
        else:
            return abs(mesh.convex_hull.volume)
    
    def compute_surface_deviation(self, baseline: trimesh.Trimesh, adjusted: trimesh.Trimesh, method: str = 'point_to_surface') -> Dict:
        if method == 'point_to_surface':
            closest_points, distances, triangle_ids = trimesh.proximity.closest_point(baseline, adjusted.vertices)
            vectors = adjusted.vertices - closest_points
            face_normals = baseline.face_normals[triangle_ids]
            signs = np.sign(np.sum(vectors * face_normals, axis=1))
            signed_distances = distances * signs
        else:
            tree = cKDTree(baseline.vertices)
            distances, indices = tree.query(adjusted.vertices, k=1)
            vectors = adjusted.vertices - baseline.vertices[indices]
            vertex_normals = baseline.vertex_normals[indices]
            signs = np.sign(np.sum(vectors * vertex_normals, axis=1))
            signed_distances = distances * signs
            closest_points = baseline.vertices[indices]
        
        self.distances = signed_distances
        self.closest_points = closest_points
        rms = np.sqrt(np.mean(signed_distances ** 2))
        adjustment_threshold = 0.01
        adjustment_mask = signed_distances < -adjustment_threshold
        vertex_areas = self._compute_vertex_areas(adjusted)
        adjustment_area = np.sum(vertex_areas[adjustment_mask])
        
        return {
            'distances': signed_distances,
            'closest_points': closest_points,
            'rms_um': rms * 1000,
            'mean_um': np.mean(signed_distances) * 1000,
            'std_um': np.std(signed_distances) * 1000,
            'max_deviation_um': np.max(signed_distances) * 1000,
            'min_deviation_um': np.min(signed_distances) * 1000,
            'median_um': np.median(signed_distances) * 1000,
            'adjustment_area_mm2': adjustment_area,
            'n_vertices_adjusted': np.sum(adjustment_mask),
            'percent_surface_adjusted': (np.sum(adjustment_mask) / len(signed_distances)) * 100
        }
    
    def _compute_vertex_areas(self, mesh: trimesh.Trimesh) -> np.ndarray:
        vertex_areas = np.zeros(len(mesh.vertices))
        face_areas = mesh.area_faces
        for i, face in enumerate(mesh.faces):
            area_per_vertex = face_areas[i] / 3
            vertex_areas[face] += area_per_vertex
        return vertex_areas
    
    def analyze(self, baseline_path: str, adjusted_path: str) -> Dict:
        self.baseline_mesh, baseline_info = self.load_mesh(baseline_path)
        self.adjusted_mesh, adjusted_info = self.load_mesh(adjusted_path)
        baseline_volume = self.calculate_volume(self.baseline_mesh)
        adjusted_volume = self.calculate_volume(self.adjusted_mesh)
        delta_volume = baseline_volume - adjusted_volume
        deviation_results = self.compute_surface_deviation(self.baseline_mesh, self.adjusted_mesh)
        return {
            'baseline': baseline_info,
            'adjusted': adjusted_info,
            'volume_analysis': {
                'baseline_volume_mm3': baseline_volume,
                'adjusted_volume_mm3': adjusted_volume,
                'delta_volume_mm3': delta_volume,
                'percent_volume_removed': (delta_volume / baseline_volume) * 100 if baseline_volume > 0 else 0
            },
            'deviation_analysis': deviation_results,
            'timestamp': datetime.now().isoformat()
        }


class AnalysisWorker(QThread):
    progress = pyqtSignal(int, str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
    
    def __init__(self, analyzer: MeshAnalyzer, baseline_path: str, adjusted_path: str):
        super().__init__()
        self.analyzer = analyzer
        self.baseline_path = baseline_path
        self.adjusted_path = adjusted_path
        
    def run(self):
        try:
            self.progress.emit(10, "Loading baseline mesh...")
            self.analyzer.baseline_mesh, _ = self.analyzer.load_mesh(self.baseline_path)
            self.progress.emit(30, "Loading adjusted mesh...")
            self.analyzer.adjusted_mesh, _ = self.analyzer.load_mesh(self.adjusted_path)
            self.progress.emit(50, "Calculating volumes...")
            self.progress.emit(70, "Computing surface deviations...")
            results = self.analyzer.analyze(self.baseline_path, self.adjusted_path)
            self.progress.emit(100, "Analysis complete!")
            self.finished.emit(results)
        except Exception as e:
            self.error.emit(str(e))


class CrownAnalyzerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.analyzer = MeshAnalyzer()
        self.baseline_path: Optional[str] = None
        self.adjusted_path: Optional[str] = None
        self.results: Optional[Dict] = None
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle("Crown Adjustment Analyzer v1.3")
        self.setMinimumSize(1400, 900)
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left panel
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        file_group = self._create_file_group()
        left_layout.addWidget(file_group)
        options_group = self._create_options_group()
        left_layout.addWidget(options_group)
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        left_layout.addWidget(self.progress_bar)
        self.progress_label = QLabel("")
        self.progress_label.setAlignment(Qt.AlignCenter)
        left_layout.addWidget(self.progress_label)
        
        self.results_tabs = QTabWidget()
        self.volume_text = QTextEdit()
        self.volume_text.setReadOnly(True)
        self.volume_text.setFont(QFont("Consolas", 10))
        self.results_tabs.addTab(self.volume_text, "Volume Analysis")
        self.deviation_text = QTextEdit()
        self.deviation_text.setReadOnly(True)
        self.deviation_text.setFont(QFont("Consolas", 10))
        self.results_tabs.addTab(self.deviation_text, "Deviation Analysis")
        self.report_text = QTextEdit()
        self.report_text.setReadOnly(True)
        self.report_text.setFont(QFont("Consolas", 9))
        self.results_tabs.addTab(self.report_text, "Full Report")
        left_layout.addWidget(self.results_tabs)
        
        export_layout = QHBoxLayout()
        self.export_csv_btn = QPushButton("Export CSV")
        self.export_csv_btn.clicked.connect(self.export_csv)
        self.export_csv_btn.setEnabled(False)
        export_layout.addWidget(self.export_csv_btn)
        self.export_report_btn = QPushButton("Export Report")
        self.export_report_btn.clicked.connect(self.export_report)
        self.export_report_btn.setEnabled(False)
        export_layout.addWidget(self.export_report_btn)
        self.export_screenshot_btn = QPushButton("Save Screenshot")
        self.export_screenshot_btn.clicked.connect(self.save_screenshot)
        self.export_screenshot_btn.setEnabled(False)
        export_layout.addWidget(self.export_screenshot_btn)
        left_layout.addLayout(export_layout)
        
        # Right panel - 3D Visualization
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        viz_group = QGroupBox("3D Visualization")
        viz_layout = QVBoxLayout(viz_group)
        
        # Row 1 controls
        viz_controls = QHBoxLayout()
        viz_controls.addWidget(QLabel("View:"))
        self.view_combo = QComboBox()
        self.view_combo.addItems(["Heat Map", "Baseline Only", "Adjusted Only", "Side by Side"])
        self.view_combo.currentIndexChanged.connect(self.update_visualization)
        viz_controls.addWidget(self.view_combo)
        viz_controls.addWidget(QLabel("Colormap:"))
        self.colormap_combo = QComboBox()
        self.colormap_combo.addItems(["RdYlBu_r", "coolwarm", "jet", "viridis", "plasma"])
        self.colormap_combo.currentIndexChanged.connect(self.update_visualization)
        viz_controls.addWidget(self.colormap_combo)
        viz_controls.addWidget(QLabel("Range (µm):"))
        self.range_spin = QSpinBox()
        self.range_spin.setRange(10, 500)
        self.range_spin.setValue(100)
        self.range_spin.valueChanged.connect(self.update_visualization)
        viz_controls.addWidget(self.range_spin)
        self.symmetric_check = QCheckBox("Symmetric")
        self.symmetric_check.setChecked(True)
        self.symmetric_check.stateChanged.connect(self.update_visualization)
        viz_controls.addWidget(self.symmetric_check)
        viz_controls.addStretch()
        viz_layout.addLayout(viz_controls)
        
        # Row 2 controls - ANATOMY VISIBILITY & BRIGHTNESS
        viz_controls2 = QHBoxLayout()
        viz_controls2.addWidget(QLabel("Opacity:"))
        self.opacity_spin = QDoubleSpinBox()
        self.opacity_spin.setRange(0.3, 1.0)
        self.opacity_spin.setValue(1.0)
        self.opacity_spin.setSingleStep(0.05)
        self.opacity_spin.setDecimals(2)
        self.opacity_spin.valueChanged.connect(self.update_visualization)
        viz_controls2.addWidget(self.opacity_spin)
        
        # Brightness slider
        viz_controls2.addWidget(QLabel("Brightness:"))
        self.brightness_slider = QSlider(Qt.Horizontal)
        self.brightness_slider.setRange(10, 100)  # 10% to 100%
        self.brightness_slider.setValue(50)  # Default 50%
        self.brightness_slider.setFixedWidth(100)
        self.brightness_slider.valueChanged.connect(self.update_visualization)
        viz_controls2.addWidget(self.brightness_slider)
        self.brightness_label = QLabel("50%")
        self.brightness_label.setFixedWidth(35)
        viz_controls2.addWidget(self.brightness_label)
        
        self.show_edges_check = QCheckBox("Show Edges")
        self.show_edges_check.setChecked(False)
        self.show_edges_check.stateChanged.connect(self.update_visualization)
        viz_controls2.addWidget(self.show_edges_check)
        
        self.silhouette_check = QCheckBox("Silhouette")
        self.silhouette_check.setChecked(False)
        self.silhouette_check.stateChanged.connect(self.update_visualization)
        viz_controls2.addWidget(self.silhouette_check)
        
        viz_controls2.addStretch()
        viz_layout.addLayout(viz_controls2)
        
        # PyVista plotter
        self.plotter_frame = QFrame()
        self.plotter_frame.setMinimumSize(600, 500)
        plotter_layout = QVBoxLayout(self.plotter_frame)
        self.plotter = QtInteractor(self.plotter_frame)
        plotter_layout.addWidget(self.plotter.interactor)
        viz_layout.addWidget(self.plotter_frame)
        right_layout.addWidget(viz_group)
        
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([500, 900])
        self._apply_styling()
        
    def _create_file_group(self) -> QGroupBox:
        group = QGroupBox("File Selection")
        layout = QGridLayout(group)
        layout.addWidget(QLabel("Baseline Crown (Pre-adjustment):"), 0, 0)
        self.baseline_label = QLabel("No file selected")
        self.baseline_label.setStyleSheet("color: gray; font-style: italic;")
        layout.addWidget(self.baseline_label, 0, 1)
        self.baseline_btn = QPushButton("Browse...")
        self.baseline_btn.clicked.connect(lambda: self.load_file('baseline'))
        layout.addWidget(self.baseline_btn, 0, 2)
        layout.addWidget(QLabel("Adjusted Crown (Post-adjustment):"), 1, 0)
        self.adjusted_label = QLabel("No file selected")
        self.adjusted_label.setStyleSheet("color: gray; font-style: italic;")
        layout.addWidget(self.adjusted_label, 1, 1)
        self.adjusted_btn = QPushButton("Browse...")
        self.adjusted_btn.clicked.connect(lambda: self.load_file('adjusted'))
        layout.addWidget(self.adjusted_btn, 1, 2)
        self.analyze_btn = QPushButton("Analyze")
        self.analyze_btn.setEnabled(False)
        self.analyze_btn.clicked.connect(self.run_analysis)
        self.analyze_btn.setStyleSheet("""
            QPushButton { background-color: #2196F3; color: white; font-weight: bold; padding: 10px; border-radius: 5px; }
            QPushButton:hover { background-color: #1976D2; }
            QPushButton:disabled { background-color: #BDBDBD; }
        """)
        layout.addWidget(self.analyze_btn, 2, 0, 1, 3)
        return group
    
    def _create_options_group(self) -> QGroupBox:
        group = QGroupBox("Analysis Options")
        layout = QGridLayout(group)
        layout.addWidget(QLabel("Distance Method:"), 0, 0)
        self.method_combo = QComboBox()
        self.method_combo.addItems(["Point-to-Surface (Recommended)", "Point-to-Point"])
        layout.addWidget(self.method_combo, 0, 1)
        layout.addWidget(QLabel("Adjustment Threshold (µm):"), 1, 0)
        self.threshold_spin = QDoubleSpinBox()
        self.threshold_spin.setRange(1, 100)
        self.threshold_spin.setValue(10)
        self.threshold_spin.setDecimals(1)
        layout.addWidget(self.threshold_spin, 1, 1)
        return group
    
    def _apply_styling(self):
        self.setStyleSheet("""
            QMainWindow { background-color: #F5F5F5; }
            QGroupBox { font-weight: bold; border: 2px solid #E0E0E0; border-radius: 5px; margin-top: 10px; padding-top: 10px; }
            QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; }
            QTextEdit { border: 1px solid #E0E0E0; border-radius: 3px; }
        """)
        
    def load_file(self, file_type: str):
        filepath, _ = QFileDialog.getOpenFileName(self, f"Select {file_type.capitalize()} Crown STL", "", "STL Files (*.stl);;All Files (*)")
        if filepath:
            if file_type == 'baseline':
                self.baseline_path = filepath
                self.baseline_label.setText(os.path.basename(filepath))
                self.baseline_label.setStyleSheet("color: green;")
            else:
                self.adjusted_path = filepath
                self.adjusted_label.setText(os.path.basename(filepath))
                self.adjusted_label.setStyleSheet("color: green;")
            if self.baseline_path and self.adjusted_path:
                self.analyze_btn.setEnabled(True)
                
    def run_analysis(self):
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.analyze_btn.setEnabled(False)
        self.worker = AnalysisWorker(self.analyzer, self.baseline_path, self.adjusted_path)
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.analysis_complete)
        self.worker.error.connect(self.analysis_error)
        self.worker.start()
        
    def update_progress(self, value: int, message: str):
        self.progress_bar.setValue(value)
        self.progress_label.setText(message)
        
    def analysis_complete(self, results: Dict):
        self.results = results
        self.progress_bar.setVisible(False)
        self.progress_label.setText("Analysis complete!")
        self.analyze_btn.setEnabled(True)
        self.export_csv_btn.setEnabled(True)
        self.export_report_btn.setEnabled(True)
        self.export_screenshot_btn.setEnabled(True)
        self.display_results()
        self.update_visualization()
        
    def analysis_error(self, error_msg: str):
        self.progress_bar.setVisible(False)
        self.progress_label.setText("")
        self.analyze_btn.setEnabled(True)
        QMessageBox.critical(self, "Analysis Error", f"An error occurred:\n{error_msg}")
        
    def display_results(self):
        if not self.results:
            return
        vol = self.results['volume_analysis']
        dev = self.results['deviation_analysis']
        
        volume_text = f"""
╔══════════════════════════════════════════════════════════════╗
║                    VOLUME ANALYSIS RESULTS                   ║
╚══════════════════════════════════════════════════════════════╝

  Baseline Crown Volume:     {vol['baseline_volume_mm3']:.4f} mm³
  Adjusted Crown Volume:     {vol['adjusted_volume_mm3']:.4f} mm³
  
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  
  DELTA VOLUME (Removed):    {vol['delta_volume_mm3']:.4f} mm³
  Percent Volume Removed:    {vol['percent_volume_removed']:.2f}%
"""
        self.volume_text.setText(volume_text)
        
        deviation_text = f"""
╔══════════════════════════════════════════════════════════════╗
║                  DEVIATION ANALYSIS RESULTS                  ║
╚══════════════════════════════════════════════════════════════╝

  RMS Deviation:             {dev['rms_um']:.2f} µm
  Mean Deviation:            {dev['mean_um']:.2f} µm
  Standard Deviation:        {dev['std_um']:.2f} µm
  Median Deviation:          {dev['median_um']:.2f} µm
  
  Maximum Deviation:         {dev['max_deviation_um']:.2f} µm
  Minimum Deviation:         {dev['min_deviation_um']:.2f} µm
  
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  
  Estimated Adjustment Area: {dev['adjustment_area_mm2']:.4f} mm²
  Surface Area Adjusted:     {dev['percent_surface_adjusted']:.2f}%
"""
        self.deviation_text.setText(deviation_text)
        
        baseline = self.results['baseline']
        adjusted = self.results['adjusted']
        report_text = f"""
================================================================================
                    CROWN ADJUSTMENT ANALYSIS REPORT
================================================================================
Generated: {self.results['timestamp']}

INPUT FILES
--------------------------------------------------------------------------------
Baseline: {baseline['filename']} ({baseline['vertices']} vertices, {baseline['faces']} faces)
Adjusted: {adjusted['filename']} ({adjusted['vertices']} vertices, {adjusted['faces']} faces)

VOLUME ANALYSIS
--------------------------------------------------------------------------------
Baseline Volume:          {vol['baseline_volume_mm3']:.6f} mm³
Adjusted Volume:          {vol['adjusted_volume_mm3']:.6f} mm³
Delta Volume (Removed):   {vol['delta_volume_mm3']:.6f} mm³
Percent Volume Removed:   {vol['percent_volume_removed']:.4f}%

SURFACE DEVIATION ANALYSIS
--------------------------------------------------------------------------------
RMS Deviation:            {dev['rms_um']:.4f} µm
Mean Deviation:           {dev['mean_um']:.4f} µm
Standard Deviation:       {dev['std_um']:.4f} µm
Max Deviation:            {dev['max_deviation_um']:.4f} µm
Min Deviation:            {dev['min_deviation_um']:.4f} µm
Adjustment Area:          {dev['adjustment_area_mm2']:.6f} mm²
================================================================================
"""
        self.report_text.setText(report_text)
        
    def update_visualization(self):
        """Update 3D visualization with adjustable brightness."""
        if not self.results or self.analyzer.adjusted_mesh is None:
            return
        self.plotter.clear()
        
        view_mode = self.view_combo.currentText()
        colormap = self.colormap_combo.currentText()
        symmetric = self.symmetric_check.isChecked()
        opacity = self.opacity_spin.value()
        show_edges = self.show_edges_check.isChecked()
        show_silhouette = self.silhouette_check.isChecked()
        
        # Get brightness value (10-100) and convert to intensity factor
        brightness = self.brightness_slider.value()
        self.brightness_label.setText(f"{brightness}%")
        
        # Setup lighting with adjustable brightness
        self._setup_even_lighting(brightness / 100.0)
        
        # Calculate ambient based on brightness (higher brightness = higher ambient)
        ambient = 0.3 + (brightness / 100.0) * 0.4  # Range: 0.3 to 0.7
        diffuse = 0.3 + (brightness / 100.0) * 0.3  # Range: 0.3 to 0.6
        
        if view_mode == "Heat Map":
            pv_mesh = pv.wrap(self.analyzer.adjusted_mesh)
            distances_um = self.results['deviation_analysis']['distances'] * 1000
            pv_mesh.point_data['Deviation (µm)'] = distances_um
            
            if symmetric:
                clim = [-self.range_spin.value(), self.range_spin.value()]
            else:
                clim = [distances_um.min(), distances_um.max()]
            
            self.plotter.add_mesh(
                pv_mesh, 
                scalars='Deviation (µm)', 
                cmap=colormap, 
                clim=clim,
                opacity=opacity, 
                show_edges=show_edges, 
                edge_color='gray', 
                line_width=0.5,
                smooth_shading=True,
                ambient=ambient,
                diffuse=diffuse,
                specular=0.0,
                show_scalar_bar=True,
                scalar_bar_args={
                    'title': 'Deviation (µm)',
                    'vertical': True,
                    'position_x': 0.88,
                    'position_y': 0.1,
                    'width': 0.06,
                    'height': 0.8,
                    'title_font_size': 14,
                    'label_font_size': 12,
                    'n_labels': 9,
                    'fmt': '%.0f',
                    'shadow': False,
                    'color': 'black'
                }
            )
            if show_silhouette:
                self.plotter.add_silhouette(pv_mesh, color='black', line_width=2)
                
        elif view_mode == "Baseline Only":
            pv_mesh = pv.wrap(self.analyzer.baseline_mesh)
            self.plotter.add_mesh(
                pv_mesh, color='lightblue', opacity=opacity, 
                show_edges=show_edges, edge_color='gray', 
                smooth_shading=True, ambient=ambient, diffuse=diffuse, specular=0.0
            )
            if show_silhouette:
                self.plotter.add_silhouette(pv_mesh, color='black', line_width=2)
                
        elif view_mode == "Adjusted Only":
            pv_mesh = pv.wrap(self.analyzer.adjusted_mesh)
            self.plotter.add_mesh(
                pv_mesh, color='lightcoral', opacity=opacity, 
                show_edges=show_edges, edge_color='gray', 
                smooth_shading=True, ambient=ambient, diffuse=diffuse, specular=0.0
            )
            if show_silhouette:
                self.plotter.add_silhouette(pv_mesh, color='black', line_width=2)
                
        elif view_mode == "Side by Side":
            bounds = self.analyzer.baseline_mesh.bounds
            offset = (bounds[1][0] - bounds[0][0]) * 0.7
            baseline_pv = pv.wrap(self.analyzer.baseline_mesh)
            baseline_pv.translate([-offset, 0, 0])
            self.plotter.add_mesh(
                baseline_pv, color='lightblue', opacity=opacity, 
                label='Baseline', show_edges=show_edges, 
                smooth_shading=True, ambient=ambient, diffuse=diffuse, specular=0.0
            )
            if show_silhouette:
                self.plotter.add_silhouette(baseline_pv, color='darkblue', line_width=2)
            adjusted_pv = pv.wrap(self.analyzer.adjusted_mesh)
            adjusted_pv.translate([offset, 0, 0])
            self.plotter.add_mesh(
                adjusted_pv, color='lightcoral', opacity=opacity, 
                label='Adjusted', show_edges=show_edges, 
                smooth_shading=True, ambient=ambient, diffuse=diffuse, specular=0.0
            )
            if show_silhouette:
                self.plotter.add_silhouette(adjusted_pv, color='darkred', line_width=2)
            self.plotter.add_legend()
        
        self.plotter.add_axes()
        self.plotter.reset_camera()
        self.plotter.update()
        
    def _setup_even_lighting(self, intensity_factor: float = 0.5):
        """Setup uniform lighting with adjustable intensity."""
        self.plotter.remove_all_lights()
        
        # Base intensity scaled by the brightness factor
        base_intensity = 0.2 + (intensity_factor * 0.5)  # Range: 0.2 to 0.7
        
        # Add lights from 6 directions for even illumination
        light_positions = [
            (1, 0, 0),    # Right
            (-1, 0, 0),   # Left
            (0, 1, 0),    # Front
            (0, -1, 0),   # Back
            (0, 0, 1),    # Top
            (0, 0, -1),   # Bottom
        ]
        
        for pos in light_positions:
            light = pv.Light(
                position=pos,
                focal_point=(0, 0, 0),
                intensity=base_intensity,
                positional=False,
                cone_angle=180,
            )
            self.plotter.add_light(light)
        
    def export_csv(self):
        if not self.results: return
        filepath, _ = QFileDialog.getSaveFileName(self, "Save CSV", f"crown_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv", "CSV Files (*.csv)")
        if filepath:
            vol, dev = self.results['volume_analysis'], self.results['deviation_analysis']
            with open(filepath, 'w') as f:
                f.write("Parameter,Value,Unit\n")
                f.write(f"Baseline File,{self.results['baseline']['filename']},\n")
                f.write(f"Adjusted File,{self.results['adjusted']['filename']},\n")
                f.write(f"Baseline Volume,{vol['baseline_volume_mm3']:.6f},mm³\n")
                f.write(f"Adjusted Volume,{vol['adjusted_volume_mm3']:.6f},mm³\n")
                f.write(f"Delta Volume,{vol['delta_volume_mm3']:.6f},mm³\n")
                f.write(f"Percent Volume Removed,{vol['percent_volume_removed']:.4f},%\n")
                f.write(f"RMS Deviation,{dev['rms_um']:.4f},µm\n")
                f.write(f"Mean Deviation,{dev['mean_um']:.4f},µm\n")
                f.write(f"Std Deviation,{dev['std_um']:.4f},µm\n")
                f.write(f"Max Deviation,{dev['max_deviation_um']:.4f},µm\n")
                f.write(f"Min Deviation,{dev['min_deviation_um']:.4f},µm\n")
                f.write(f"Adjustment Area,{dev['adjustment_area_mm2']:.6f},mm²\n")
                f.write(f"Percent Surface Adjusted,{dev['percent_surface_adjusted']:.4f},%\n")
            QMessageBox.information(self, "Export Complete", f"CSV saved to:\n{filepath}")
            
    def export_report(self):
        if not self.results: return
        filepath, _ = QFileDialog.getSaveFileName(self, "Save Report", f"crown_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt", "Text Files (*.txt)")
        if filepath:
            with open(filepath, 'w') as f:
                f.write(self.report_text.toPlainText())
            QMessageBox.information(self, "Export Complete", f"Report saved to:\n{filepath}")
            
    def save_screenshot(self):
        if not self.results: return
        filepath, _ = QFileDialog.getSaveFileName(self, "Save Screenshot", f"crown_heatmap_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png", "PNG Files (*.png)")
        if filepath:
            self.plotter.screenshot(filepath)
            QMessageBox.information(self, "Screenshot Saved", f"Screenshot saved to:\n{filepath}")
            
    def closeEvent(self, event):
        self.plotter.close()
        event.accept()


def main():
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication(sys.argv)
    app.setApplicationName("Crown Adjustment Analyzer")
    window = CrownAnalyzerGUI()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
