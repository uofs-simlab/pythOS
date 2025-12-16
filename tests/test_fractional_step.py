import pytest
import numpy as np
import sys
sys.path.insert(1,sys.path[0]+"/../")
import fractional_step as fs
fs.alphas_repo([], 0, 2)
#'Godunov', 'Godunov-3', 'Godunov-5', 'Godunov-N'
#'SM2', 'Strang', 'Strang-3', 'StrangBCA-3', 'StrangBAC-3', 'StrangCBA-3', 'StrangCAB-3', 'StrangACB-3', 'Strang-3split', 'StrangACB-3split', 
#'Best22', 'AKOpt22', 'Milne22_Complex_i','Strang-5','Strang-N'
#'R3', 'C3','S3','AKS3','AKS3C', 'AKS3P','PP3_4A-3', 'OS43(7)'
#'Y4', 'Y4-3', 'M4', 'B4',  'C4', 'C4-3', 'Yoshida-3', 'Yoshida-4',  'Yoshida-5'
class TestExact:
    @pytest.mark.parametrize("method", ["Godunov", "Godunov-N"])
    def test_order_1(self,method):
        A = 2
        B = 1
        def f1(t,y):
            return np.ones(y.size)*A
        def f2(t,y):
            return np.array([B,0])
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A + np.array([B,0]))*tf + y0)
        assert (fs.fractional_step([f1,f2],dt,y0,t0,tf,method,{(0,):"ADAPTIVE"}) == yf)
    @pytest.mark.parametrize("method", ["Godunov-3", "Godunov-N"])
    def test_order_1_N3(self,method):
        A = 2
        B = 1
        C = -0.5
        def f1(t,y):
            return np.ones(y.size)*A
        def f2(t,y):
            return np.array([B,0])
        def f3(t, y):
            return np.array([0, C])
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A + np.array([B,C]))*tf + y0)
        assert (fs.fractional_step([f1,f2, f3],dt,y0,t0,tf,method,{(0,):"ADAPTIVE"}) == yf)
    
    @pytest.mark.parametrize("method", ['Strang-3', 'StrangBCA-3', 'StrangBAC-3', 'StrangCBA-3', 'StrangCAB-3', 'StrangACB-3', 'Strang-3split', 'StrangACB-3split','Strang-N'])
    def test_order_2_N3(self,method):
        A1 = np.array([1,0])
        B1 = -0.5
        def f1(t,y):
            return (A1 * t + B1)
        A2 = np.array([1,5])
        B2 = 1
        def f2(t,y):
            return (A2 * t + B2)
        A3 = 0.5
        B3 = np.array([4,1])
        def f3(t,y):
            return (A3 * t + B3)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A1+A2+A3)*tf**2/2 + (B1+B2+B3)*tf + y0)
        assert (fs.fractional_step([f1,f2,f3],dt,y0,t0,tf,method,{(0,):"ADAPTIVE"}) == yf)
    @pytest.mark.parametrize("method", ['SM2', 'Strang','Best22', 'AKOpt22', 'Milne22_Complex_i','Strang-N'])
    def test_order_2(self,method):
        A1 = np.array([1,0])
        B1 = -0.5
        def f1(t,y):
            return (A1 * t + B1)
        A2 = np.array([1,5])
        B2 = 1
        def f2(t,y):
            return (A2 * t + B2)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A1+A2)*tf**2/2 + (B1+B2)*tf + y0)
        assert (fs.fractional_step([f1,f2],dt,y0,t0,tf,method,{(0,):"ADAPTIVE"},ivp_methods={1: ("Dormand-Prince", 1e-8, 1e-10),2: ("Dormand-Prince", 1e-8, 1e-10)}) == yf)

    @pytest.mark.parametrize("method", ['R3', 'C3','S3','AKS3','AKS3C', 'AKS3P','PP3_4A-3', 'OS43(7)'])
    def test_order_3(self,method):
        A1 = np.array([1,1])
        B1 = -0.5
        C1 = np.array([2,3])
        def f1(t,y):
            return (A1 * t*t + B1*t+C1)
        A2 = np.array([1,4])
        B2 = 1.5
        C2 = np.array([1,0])
        def f2(t,y):
            return (A2 * t*t + B2*t+C2)

        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A1+A2)*tf**3/3 + (B1+B2)*tf**2/2 + (C1+C2)*tf + y0)
        assert (fs.fractional_step([f1,f2],dt,y0,t0,tf,method,{(0,):"ADAPTIVE"}, ivp_methods={1: ("Dormand-Prince", 1e-8, 1e-10),2: ("Dormand-Prince", 1e-8, 1e-10)}) == yf)

    @pytest.mark.parametrize("method", ['Y4', 'M4', 'B4',  'C4'])
    def test_order_4(self,method):
        A1 = np.array([1,0])
        B1 = -0.5
        C1 = np.array([2,3])
        D1 = np.array([-1,2])
        def f1(t,y):
            return (A1 * t**3 + B1*t**2+C1*t+D1)
        A2 = np.array([1,5])
        B2 = 1
        C2 = np.array([1,0])
        D2 = np.array([2,4])
        def f2(t,y):
            return (A2 * t**3 + B2*t**2+C2*t+D2)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A1+A2)*tf**4/4 + (B1+B2)*tf**3/3 + (C1+C2)*tf**2/2 + (D1+D2)*tf + y0)
        assert (fs.fractional_step([f1,f2],dt,y0,t0,tf,method,{(0,):"ADAPTIVE"}, ivp_methods={1: ("Dormand-Prince", 1e-8, 1e-10),2: ("Dormand-Prince", 1e-8, 1e-10)}) == yf)

    @pytest.mark.parametrize("method", ['Y4-3','C4-3', 'Yoshida-3'])
    def test_order_4_N3(self,method):
        A1 = np.array([1,0])
        B1 = -0.5
        C1 = np.array([2,3])
        D1 = np.array([-1,2])
        def f1(t,y):
            return (A1 * t**3 + B1*t**2+C1*t+D1)
        A2 = np.array([1,5])
        B2 = 1
        C2 = np.array([1,0])
        D2 = np.array([2,4])
        def f2(t,y):
            return (A2 * t**3 + B2*t**2+C2*t+D2)
        A3 = np.array([5,2])
        B3 = 1.5
        C3 = np.array([1,4])
        D3 = np.array([2.5,-1])
        def f3(t,y):
            return (A3 * t**3 + B3*t**2+C3*t+D3)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A1+A2+A3)*tf**4/4 + (B1+B2+B3)*tf**3/3 + (C1+C2+C3)*tf**2/2 + (D1+D2+D3)*tf + y0)
        assert (fs.fractional_step([f1,f2,f3],dt,y0,t0,tf,method,{(0,):"ADAPTIVE"}, ivp_methods={1: ("Dormand-Prince", 1e-8, 1e-10),2: ("Dormand-Prince", 1e-8, 1e-10),3:("Dormand-Prince", 1e-8, 1e-10)}) == yf)

class TestOrder:
    # Basic Test problem
    def f1(t,y):
        x = y[0]
        y = y[1]
        return np.array([(-y*t)/t, (x*t)/t])
    def f1_a(t,y):
        x = y[0]
        y = y[1]
        return np.array([(-y*t)/t, 0])
    def f1_b(t,y):
        x = y[0]
        y = y[1]
        return np.array([0, (x*t)/t])
    def f2(t,y):
        x = y[0]
        y = y[1]
        return np.array([(+x)/t, (-y)/t])
    
    def f2_a(t, y):
        x = y[0]
        y = y[1]
        return np.array([(+x)/t, 0])
    def f2_b(t,y):
        x = y[0]
        y = y[1]
        return np.array([0, (-y)/t])

    y0 = np.array([np.pi,-np.pi])
    t0 = np.pi
    tf = 2*np.pi
    yf = np.array([-np.pi,1/2+np.pi])

# ['R3', 'C3','S3','AKS3','AKS3C', 'AKS3P','PP3_4A-3', 'OS43(7)']
# ['Y4', 'M4', 'B4',  'C4']
    @pytest.mark.parametrize("method,subint,order,dt0,f_list", [('Godunov', {(1,): "FE", (2,): "BE"}, 1,tf/8,[f1,f2]),
                                                                ('Godunov-N', {(1,): "BE", (2,): "FE"}, 1, tf/8, [f1, f2]),
                                                                ('Godunov-N', {(1,): "BE", (3,): "BE"}, 1, tf/8, [f1_a, f2_a, f1_b, f2_b]),
                                                                ('SM2',  {(1,): "Heun", (2,): "GL2"}, 2, tf/8, [f1, f2]), 
                                                                 ('Strang', {(1,): "Heun", (2,): "GL2"}, 2, tf/8, [f1, f2]), 
                                                                 ('Best22',  {(1,): "Heun", (2,): "GL2"}, 2, tf/8, [f1, f2]), 
                                                                 ('AKOpt22', {(1,): "Heun", (2,): "GL2"}, 2, tf/8, [f1, f2]), 
                                                                 ('Milne22_Complex_i', {(1,): "Heun", (2,): "GL2"}, 2, tf/8, [f1, f2]), 
                                                                 ('Strang-N', {(1,): "Heun", (2,): "GL2"}, 2, tf/8, [f1, f2]),
                                                                 ('Strang-3', {(1,): "Heun", (2,): "Heun", (3,): "Heun"}, 2, tf/8, [f1_a, f1_b, f2]),
                                                                 ('Strang-N', {(1,): "Heun", (2,): "GL2", (3,): "Heun", (4,): "Heun"}, 2, tf/8, [f1_a, f1_b, f2_a, f2_b]),
                                                                 ('R3', {(1,): "RK3", (2,): "SD2O3"}, 3, tf/8, [f1, f2]),
                                                                 ('C3', {(1,): "RK3", (2,): "SD2O3"}, 3, tf/8, [f1, f2]),
                                                                 ('S3', {(1,): "RK3", (2,): "SD2O3"}, 3, tf/8, [f1, f2]),
                                                                 ('AKS3C', {(1,): "RK3", (2,): "SD2O3"}, 3, tf/8, [f1, f2]),
                                                                 ('AKS3P', {(1,): "RK3", (2,): "SD2O3"}, 3, tf/8, [f1, f2]),
                                                                 ('Y4', {(1,): "RK4", (2,): "SD3O4"}, 4, tf/8, [f1, f2]),
                                                                 ('M4', {(1,): "RK4", (2,): "SD3O4"}, 4, tf/8, [f1, f2]),
                                                                 ('C4', {(1,): "RK4", (2,): "SD3O4"}, 4, tf/8, [f1, f2])])
    def test_order(self,method,subint,order,dt0,f_list):
        for dt in [dt0,dt0/2,dt0/4]:
            assert np.log2(np.linalg.norm(fs.fractional_step(f_list,dt,TestOrder.y0,TestOrder.t0,TestOrder.tf,method,subint)-TestOrder.yf)/(np.linalg.norm(fs.fractional_step(f_list,dt/2,TestOrder.y0,TestOrder.t0,TestOrder.tf,method,subint)-TestOrder.yf))) == pytest.approx(order, abs=0.4, rel=0.1)

