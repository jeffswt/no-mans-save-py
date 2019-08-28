
import json
import re
import random


class vec:
    def __init__(self, *args):
        if len(args) == 0:
            self.x, self.y, self.z = 0.0, 0.0, 0.0
        elif len(args) == 1:
            if type(args[0]) in (int, float):
                fv = float(args[0])
                self.x, self.y, self.z = fv, fv, fv
            elif type(args[0]) in (tuple, list):
                ar = (float(args[0][0]), float(args[0][1]), float(args[0][2]))
                self.x, self.y, self.z = ar
            elif type(args[0]) in (vec,):
                self.x, self.y, self.z = args[0].x, args[0].y, args[0].z
            else:
                raise TypeError('unknown type')
        elif len(args) == 3:
            ar = (float(args[0]), float(args[1]), float(args[2]))
            self.x, self.y, self.z = ar
        else:
            raise AttributeError('vec() takes 0, 1, 3 arguments')

    def __repr__(self):
        return 'vec(%.3f, %.3f, %.3f)' % (self.x, self.y, self.z)

    def length(self):
        return (self.x ** 2 + self.y ** 2 + self.z ** 2) ** 0.5

    def __add__(self, value):
        value = vec(value)
        return vec(self.x + value.x, self.y + value.y, self.z + value.z)

    def __radd__(self, value):
        return self.__add__(value)

    def __sub__(self, value):
        value = vec(value)
        return vec(self.x - value.x, self.y - value.y, self.z - value.z)

    def __rsub__(self, value):
        value = vec(value)
        return vec(value.x - self.x, value.y - self.y, value.z - self.z)

    @staticmethod
    def dot(a, b):
        if type(a) != vec:
            a, b = b, a
        if type(a) != vec:
            raise TypeError('must have vector')
        if type(b) in (int, float):
            return vec(a.x * b, a.y * b, a.z * b)
        elif type(b) in (vec,):
            return a.x * b.x + a.y * b.y + a.z * b.z
        raise TypeError('must be either vector or real number')

    @staticmethod
    def cross(a, b):
        if type(a) != vec or type(b) != vec:
            raise TypeError('must be two vectors')
        return vec(a.y * b.z - a.z * b.y,
                   a.z * b.x - a.x * b.z,
                   a.x * b.y - a.y * b.x)

    def __mul__(self, value):
        return vec.dot(self, value)

    def __rmul__(self, value):
        return vec.dot(self, value)

    def __truediv__(self, value):
        if type(value) not in (int, float):
            raise TypeError('divisor must be real number')
        value = float(value)
        return vec(self.x / value, self.y / value, self.z / value)

    @staticmethod
    def norm(a):
        if type(a) != vec:
            raise TypeError('must be a vector')
        return a / a.length()

    def normalize(self):
        l = self.length()
        self.x, self.y, self.z = self.x / l, self.y / l, self.z / l
        return vec(self)
    pass


def rps(arr):
    n = len(arr)
    res = 0
    for i in range(0, n):
        for j in range(i + 1, n):
            if arr[i] > arr[j]:
                res += 1
    return res


def det(arr):
    n = len(arr)
    if n == 0:
        return 0
    for row in arr:
        if len(row) != n:
            raise ValueError('must be valid determinant')

    def itert(n, arr, depth, rows):
        if depth == n:
            res = (-1) ** rps(rows)
            for i in range(0, n):
                res *= arr[i][rows[i]]
            return res
        res = 0
        for i in range(0, n):
            if i in rows:
                continue
            res += itert(n, arr, depth + 1, rows + [i])
        return res
    return itert(n, arr, 0, [])


def matinv(m):
    n = len(m)
    adet = det(m)
    if float(adet) == 0.0:
        return None
    res = [[0.0 for i in range(0, n)] for j in range(0, n)]
    for i in range(0, n):
        for j in range(0, n):
            rm = [[m[a if a < i else a + 1][b if b < j else b + 1]
                   for b in range(0, n - 1)] for a in range(0, n - 1)]
            res[j][i] = (- 1) ** (i + j) * det(rm) / adet
    return res


def mattranspose(m):
    n = len(m)
    m_ = [[m[j][i] for j in range(0, n)] for i in range(0, n)]
    return m_


def solve(A, X):
    res = []
    Ainv = matinv(A)
    if Ainv is None:
        return None
    n = len(A)
    for i in range(0, n):
        t = 0.0
        for j in range(0, n):
            t += Ainv[i][j] * X[j]
        res.append(t)
    return res


class CoordSystem:
    def __init__(self, *args):
        if len(args) == 0:
            self.O = vec(0.0, 0.0, 0.0)
            self.X = vec(1.0, 0.0, 0.0)
            self.Y = vec(0.0, 1.0, 0.0)
            self.Z = vec(0.0, 0.0, 1.0)
        elif len(args) == 4:
            self.O = vec(args[0])
            self.X = vec.norm(vec(args[1]))
            self.Y = vec.norm(vec(args[2]))
            self.Z = vec.norm(vec(args[3]))
        else:
            raise AttributeError('CoordSystem() takes 0, 4 arguments')

    def apply(self, value):
        # map relative coords to actual coords
        v = vec(value)
        return self.O + v.x * self.X + v.y * self.Y + v.z * self.Z

    def solve(self, value):
        # map actual coords to relative coords
        val = value - self.O
        val = [val.x, val.y, val.z]
        mat = [[self.X.x, self.Y.x, self.Z.x],
               [self.X.y, self.Y.y, self.Z.y],
               [self.X.z, self.Y.z, self.Z.z]]
        return vec(solve(mat, val))
    pass


def move_base(save_filename, target_filename, base_patt, swap_name):
    """move_base(sf, tf, bn, sid) -- Move base computer to a target
    @param sf: Filename to the original save file (.hg)
    @param tf: Filename to the target save file (.hg)
    @param bn: A RegEx pattern describing the name of the base
    @param sid: The unique object (ID) to swap the base computer with"""
    swap_id = {
        'Signal Booster': '^BUILDSIGNAL',
        'Save Beacon': '^BUILDBEACON',
        'Save Point': '^BUILDSAVE',
    }.get(swap_name, swap_name)
    # Load save file
    save_file = open(save_filename, 'rb')
    save_data_raw = json.loads(save_file.read()[:-1].decode('utf-8'))
    save_file.close()
    # Parse json map
    json_map_f = open('jsonmap.txt', 'r', encoding='utf-8')
    json_map = {}
    json_map_r = {}
    for i in json_map_f.read().strip().split('\n'):
        _ = i.split('\t')
        if len(_) == 1:
            _.append(_[0])
        json_map[_[0]] = _[1]
        json_map_r[_[1]] = _[0]

    # Convert json
    def conv(s, mp):
        if isinstance(s, str):
            return s
        elif isinstance(s, list):
            return list(conv(i, mp) for i in s)
        elif isinstance(s, dict):
            return dict((mp.get(i, i), conv(s[i], mp)) for i in s)
        return s
    save_data = conv(save_data_raw, json_map)
    # Locate base
    hbase = []
    for i in save_data['PlayerStateData']['PersistentPlayerBases']:
        if re.findall(base_patt, i['Name']):
            hbase.append(i)
    if len(hbase) != 1:
        raise AttributeError('no base found under search criteria'
                             if len(hbase) == 0 else
                             'multiple bases found, narrow keywords')
    hbase = hbase[0]
    # Locate base computer and save beacon
    bcomp = []
    sboos = []
    for i in hbase['Objects']:
        oid = i['ObjectID']
        if oid == '^BASE_FLAG':
            bcomp.append(i)
        elif oid == swap_id:
            sboos.append(i)
    if len(bcomp) != 1:
        raise KeyError('you need exactly 1 base computer')
    if len(sboos) != 1:
        raise KeyError('you need exactly 1 target flag')
    bcomp = bcomp[0]
    sboos = sboos[0]
    # Extract position
    bpos = vec(hbase['Position'])
    bfrw = vec(hbase['Forward'])
    bpos_n, bfrw_n = vec.norm(bpos), vec.norm(bfrw)
    # Generate old plane
    axis_z = bfrw_n
    axis_y = vec.norm(bpos_n - vec.dot(bpos_n, bfrw_n) * bfrw_n)  # Schmidt
    axis_x = vec.norm(vec.cross(axis_y, axis_z))
    c_plane = CoordSystem(bpos, axis_x, axis_y, axis_z)
    # Generate new plane
    axis_o = c_plane.apply(vec(sboos['Position']))
    axis_y = vec.norm(axis_o)
    # nx = random.random()
    # ny = random.random()
    nx = bfrw.x
    ny = bfrw.y
    nz = - (nx * axis_y.x + ny * axis_y.y) / axis_y.z
    axis_z = vec.norm(vec(nx, ny, nz))
    axis_x = vec.norm(vec.cross(axis_y, axis_z))
    n_plane = CoordSystem(axis_o, axis_x, axis_y, axis_z)
    # Apply positions
    hbase['Position'] = [axis_o.x, axis_o.y, axis_o.z]
    hbase['Forward'] = [axis_z.x, axis_z.y, axis_z.z]
    for obj in hbase['Objects']:
        vp = vec(obj['Position'])
        vu = vec(obj['Up'])
        va = vec(obj['At'])
        ratio = 100.0
        v1 = vp
        v2 = v1 + ratio * vu
        v3 = v1 + ratio * va
        v1 = c_plane.apply(v1)
        v2 = c_plane.apply(v2)
        v3 = c_plane.apply(v3)
        v1 = n_plane.solve(v1)
        v2 = n_plane.solve(v2)
        v3 = n_plane.solve(v3)
        # ovp, ovu, ova = vec(vp), vec(vu), vec(va)
        vp = v1
        vu = (v2 - v1) / ratio
        va = (v3 - v1) / ratio
        obj['Position'] = [vp.x, vp.y, vp.z]
        obj['Up'] = [vu.x, vu.y, vu.z]
        obj['At'] = [va.x, va.y, va.z]
    # Swap base and flag
    bcomp['Position'], sboos['Position'] = sboos['Position'], bcomp['Position']
    bcomp['Up'], sboos['Up'] = sboos['Up'], bcomp['Up']
    bcomp['At'], sboos['At'] = sboos['At'], bcomp['At']
    # Save file
    save_data_raw = conv(save_data, json_map_r)
    t = open(target_filename, 'wb')
    t.write(json.dumps(save_data_raw, indent=None,
            separators=(',', ':')).encode('utf-8') + b'\x00')
    t.close()
    return

if __name__ == '__main__':
    fn = 'save.hg'
    cn = r'Colony$'
    move_base(fn, fn, cn, 'Signal Booster')
