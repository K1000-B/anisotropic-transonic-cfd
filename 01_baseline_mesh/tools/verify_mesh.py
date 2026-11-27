"""Inspect a Gmsh baseline mesh and export it natively to SU2; never run CFD."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
EXPECTED = ROOT / '.PROJECT_CFDvenv'
if Path(sys.prefix) != EXPECTED or Path(sys.executable).parent != EXPECTED / 'bin':
    raise SystemExit('Use the designated project venv interpreter.')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    # Reuse the binding bundled with the installed Gmsh; no package installation.
    binding = Path(shutil.which('gmsh')).resolve().parents[1] / 'lib' / 'gmsh.py'
    if not binding.is_file():
        raise SystemExit('Installed Gmsh binding not found; review local installation.')
    spec = importlib.util.spec_from_file_location('gmsh', binding)
    gmsh = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gmsh)
    gmsh.initialize()
    try:
        gmsh.open(str(args.mesh))
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        xy = np.asarray(coordinates).reshape(-1, 3)[:, :2]
        index = {int(tag): i for i, tag in enumerate(node_tags)}
        types, tags, nodes = gmsh.model.mesh.getElements(2)
        assert list(types) == [2], 'Only linear triangles are permitted.'
        tri = np.array([index[int(t)] for t in nodes[0]]).reshape(-1, 3)
        p = xy[tri]
        vectors = np.roll(p, -1, axis=1) - p
        lengths = np.linalg.norm(vectors, axis=2)
        cross = vectors[:, 0, 0] * (-vectors[:, 2, 1]) - vectors[:, 0, 1] * (-vectors[:, 2, 0])
        area = cross / 2
        assert np.all(area > 0), 'Inverted or degenerate triangle.'
        assert len(np.unique(np.round(xy, 12), axis=0)) == len(xy), 'Coincident nodes at 1e-12 length rounding.'
        used = np.unique(tri)
        assert len(used) == len(xy), 'Unused exported nodes.'
        edges = np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
        unique_edges, inverse, counts = np.unique(edges, axis=0, return_inverse=True, return_counts=True)
        assert np.all((counts == 1) | (counts == 2)), 'Non-manifold edge.'
        parent = np.arange(len(xy))
        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        for a, b in unique_edges:
            parent[find(a)] = find(b)
        components = len({find(i) for i in used})
        assert components == 1, 'Disconnected mesh.'
        assert len(xy)-len(unique_edges)+len(tri) == 0, 'Expected one connected domain with one hole.'
        boundary = set(map(tuple, unique_edges[counts == 1]))
        markers = {}
        marked = []
        airfoil_nodes = set()
        for dim, group in gmsh.model.getPhysicalGroups(1):
            name = gmsh.model.getPhysicalName(dim, group)
            lines = []
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, group):
                et, _, en = gmsh.model.mesh.getElements(1, int(entity))
                assert list(et) == [1]
                lines.extend(np.array([index[int(t)] for t in en[0]]).reshape(-1, 2).tolist())
            markers[name] = len(lines)
            marked.extend(tuple(sorted(e)) for e in lines)
            if name == 'airfoil':
                airfoil_nodes.update(np.asarray(lines).ravel().tolist())
        assert set(markers) == {'airfoil', 'farfield', 'outlet'}
        assert len(marked) == len(set(marked)) and set(marked) == boundary
        for name in ('airfoil', 'outer'):
            selected = [e for e in boundary if (e[0] in airfoil_nodes) == (name == 'airfoil')]
            degree = np.bincount(np.array(selected).ravel(), minlength=len(xy))
            assert np.all(degree[degree > 0] == 2), 'Boundary is not a closed conforming loop.'
        # Compare to the analytical standard NACA0012 in this default unit-chord case.
        af = xy[sorted(airfoil_nodes)]
        assert np.isclose(af[:, 0].min(), 0) and np.isclose(af[:, 0].max(), 1), 'This report expects default geometry.'
        surface = af[af[:, 0] < 1-1e-12]
        x = np.clip(surface[:, 0], 0, 1)
        yt = 0.6*(0.2969*np.sqrt(x)-0.126*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
        geometry_error = np.max(np.abs(np.abs(surface[:, 1])-yt))
        assert geometry_error < 1e-5, 'Airfoil geometry error exceeds 1e-5 chord.'
        # Sample the analytical surface and measure its distance to mesh chords.
        angle = np.linspace(0, np.pi, 20001)
        exact_x = (1-np.cos(angle))/2
        exact_y = 0.6*(0.2969*np.sqrt(exact_x)-0.126*exact_x-0.3516*exact_x**2+0.2843*exact_x**3-0.1015*exact_x**4)
        positive = af[af[:,1] >= 0]
        positive = positive[np.argsort(positive[:,0])]
        positive = positive[np.r_[True,np.diff(positive[:,0]) > 1e-12]]
        interval = np.clip(np.searchsorted(positive[:,0],exact_x)-1,0,len(positive)-2)
        a,b = positive[interval],positive[interval+1]
        d = b-a
        point = np.column_stack((exact_x,exact_y))
        t = np.clip(np.sum((point-a)*d,axis=1)/np.sum(d*d,axis=1),0,1)
        faceting_error = float(np.max(np.linalg.norm(point-a-t[:,None]*d,axis=1)))
        q = 4*np.sqrt(3)*area/np.sum(lengths**2, axis=1)
        radius = 2*area/lengths.sum(axis=1)
        aspect = lengths.max(axis=1)/(2*np.sqrt(3)*radius)
        angles = []
        for i in range(3):
            a, b, c = lengths[:, i], lengths[:, (i+1)%3], lengths[:, (i+2)%3]
            angles.append(np.degrees(np.arccos(np.clip((b*b+c*c-a*a)/(2*b*c), -1, 1))))
        adjacency = [[] for _ in unique_edges]
        for cell, edge_id in enumerate(inverse):
            adjacency[edge_id].append(cell % len(tri))
        jumps = []
        for adj in adjacency:
            if len(adj) == 2:
                a,b = adj
                jumps.append(float(np.sqrt(max(area[a],area[b])/min(area[a],area[b]))))
        def stats(values):
            return dict(zip(('min','p01','median','p95','p99','max'), map(float,np.percentile(values,[0,1,50,95,99,100]))))
        su2 = args.output/'baseline.su2'
        gmsh.write(str(su2))
        text = su2.read_text()
        exported_markers = [line.split('=',1)[1].strip() for line in text.splitlines() if line.startswith('MARKER_TAG=')]
        assert set(exported_markers) == set(markers), 'Unexpected SU2 marker.'
        assert 'NDIME= 2' in text
        lines = text.splitlines()
        ei = next(i for i,s in enumerate(lines) if s.startswith('NELEM='))
        ni = next(i for i,s in enumerate(lines) if s.startswith('NPOIN='))
        assert int(lines[ei].split('=')[1]) == len(tri)
        assert int(lines[ni].split('=')[1]) == len(xy)
        exported_tri = np.array([list(map(int,line.split()[:4])) for line in lines[ei+1:ei+1+len(tri)]])
        exported_xy = np.array([list(map(float,line.split()[:2])) for line in lines[ni+1:ni+1+len(xy)]])
        assert np.all(exported_tri[:,0] == 5), 'SU2 volume must use triangle type 5.'
        # Native export renumbers nodes; compare through coordinate identities.
        coordinate_index = {tuple(np.round(p,10)):i for i,p in enumerate(xy)}
        mapping = np.array([coordinate_index[tuple(np.round(p,10))] for p in exported_xy])
        assert len(np.unique(mapping)) == len(xy)
        assert np.allclose(exported_xy,xy[mapping],rtol=0,atol=1e-10)
        mapped = mapping[exported_tri[:,1:]]
        assert sorted(map(tuple,np.sort(mapped,axis=1))) == sorted(map(tuple,np.sort(tri,axis=1)))
        exported_boundary = []
        for i,line in enumerate(lines):
            if line.startswith('MARKER_ELEMS='):
                n = int(line.split('=')[1])
                assert all(s.split()[0]=='3' for s in lines[i+1:i+1+n])
                exported_boundary.extend(tuple(sorted(mapping[list(map(int,s.split()[1:3]))])) for s in lines[i+1:i+1+n])
        assert len(exported_boundary) == len(boundary) and set(exported_boundary) == boundary
        assert q.min() > 0.5 and min(map(np.min,angles)) > 20, 'Baseline quality acceptance failed.'
        report = {'nodes':len(xy), 'triangles':len(tri), 'markers':markers,
                  'connected_components':components, 'euler_characteristic':0,
                  'nonpositive_areas':0, 'duplicate_nodes_rounded_1e12':0,
                  'edge_length':stats(lengths.ravel()), 'area':stats(area),
                  'mean_ratio_quality':stats(q), 'minimum_angle_degrees':stats(np.min(angles,axis=0)),
                  'aspect_ratio':stats(aspect), 'neighbor_sqrt_area_ratio':stats(jumps),
                  'airfoil_boundary_max_vertical_error_over_chord':float(geometry_error),
                  'su2_markers':exported_markers, 'analytic_to_upper_polyline_sampled_distance_over_chord':faceting_error,
                  'files_sha256':{args.mesh.name:hashlib.sha256(args.mesh.read_bytes()).hexdigest(), su2.name:hashlib.sha256(su2.read_bytes()).hexdigest()},
                  'gmsh_version':gmsh.option.getString('General.Version')}
        (args.output/'quality.json').write_text(json.dumps(report,indent=2)+'\n')
        segments = xy[unique_edges]
        for name, bounds in [('overall',(-21,31,-21,21)),('airfoil',(-0.15,1.2,-0.3,0.3)),('leading_edge',(-0.04,0.12,-0.07,0.07)),('trailing_edge',(0.94,1.08,-0.045,0.045)),('downstream',(0.8,5,-1,1))]:
            fig,ax=plt.subplots(figsize=(10,6), layout='constrained')
            ax.add_collection(LineCollection(segments, colors='#314b62',linewidths=0.3))
            ax.set(xlim=bounds[:2],ylim=bounds[2:],xlabel='x/c',ylabel='y/c',title='NACA 0012 baseline mesh — '+name.replace('_',' '),aspect='equal')
            fig.savefig(args.output/(name+'.png'),dpi=180)
            plt.close(fig)
        print(json.dumps(report,indent=2))
    finally:
        gmsh.finalize()

if __name__ == '__main__':
    main()
