"""
Test suite for Crown Adjustment Analyzer
Run with: pytest tests/ -v
"""

import numpy as np
import pytest
import tempfile
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from crown_analyzer import MeshAnalyzer


def create_cube_stl(filepath, size=10.0):
    """Create a cube STL for testing."""
    import trimesh
    vertices = np.array([
        [0, 0, 0], [size, 0, 0], [size, size, 0], [0, size, 0],
        [0, 0, size], [size, 0, size], [size, size, size], [0, size, size]
    ])
    faces = np.array([
        [0, 1, 2], [0, 2, 3], [4, 6, 5], [4, 7, 6],
        [0, 4, 5], [0, 5, 1], [2, 6, 7], [2, 7, 3],
        [0, 3, 7], [0, 7, 4], [1, 5, 6], [1, 6, 2]
    ])
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    mesh.export(filepath)


def create_sphere_stl(filepath, radius=5.0):
    """Create a sphere STL for testing."""
    import trimesh
    mesh = trimesh.creation.icosphere(subdivisions=3, radius=radius)
    mesh.export(filepath)


class TestMeshAnalyzer:
    
    @pytest.fixture
    def analyzer(self):
        return MeshAnalyzer()
    
    def test_load_mesh(self, analyzer):
        """Test mesh loading."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "cube.stl")
            create_cube_stl(filepath, size=10.0)
            mesh, info = analyzer.load_mesh(filepath)
            assert info['vertices'] == 8
            assert info['faces'] == 12
    
    def test_volume_cube(self, analyzer):
        """Test volume calculation for known cube."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "cube.stl")
            create_cube_stl(filepath, size=10.0)
            mesh, _ = analyzer.load_mesh(filepath)
            volume = analyzer.calculate_volume(mesh)
            assert abs(volume - 1000.0) < 10  # 10mm cube = 1000mm³
    
    def test_identical_meshes(self, analyzer):
        """Test that identical meshes show zero deviation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "sphere.stl")
            create_sphere_stl(filepath, radius=5.0)
            mesh, _ = analyzer.load_mesh(filepath)
            results = analyzer.compute_surface_deviation(mesh, mesh)
            assert results['rms_um'] < 0.01
    
    def test_volume_removed(self, analyzer):
        """Test detection of volume reduction."""
        with tempfile.TemporaryDirectory() as tmpdir:
            baseline = os.path.join(tmpdir, "baseline.stl")
            adjusted = os.path.join(tmpdir, "adjusted.stl")
            create_sphere_stl(baseline, radius=5.0)
            create_sphere_stl(adjusted, radius=4.5)
            results = analyzer.analyze(baseline, adjusted)
            assert results['volume_analysis']['delta_volume_mm3'] > 0
            assert results['volume_analysis']['percent_volume_removed'] > 0
    
    def test_analysis_output(self, analyzer):
        """Test that analysis returns all required metrics."""
        with tempfile.TemporaryDirectory() as tmpdir:
            baseline = os.path.join(tmpdir, "baseline.stl")
            adjusted = os.path.join(tmpdir, "adjusted.stl")
            create_sphere_stl(baseline, radius=5.0)
            create_sphere_stl(adjusted, radius=4.8)
            results = analyzer.analyze(baseline, adjusted)
            
            assert 'volume_analysis' in results
            assert 'deviation_analysis' in results
            assert 'rms_um' in results['deviation_analysis']
            assert 'delta_volume_mm3' in results['volume_analysis']


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
