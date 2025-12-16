import numpy as np
import pytest
import sys
sys.path.insert(1,sys.path[0]+"/../")
from butcher_tableau import Tableau, EmbeddedTableau
from additive_rk import ark_solve

tableau_1_e = Tableau(
    np.array([0,1]),
    np.array([[0,0],[1,0]]),
    np.array([1,0]))
tableau_1_i = Tableau(
    np.array([0,1]),
    np.array([[0,0],[0,1]]),
    np.array([0,1]))
tableau_3_e = Tableau(
    np.array([0, 1767732205903/2027836641118, 3/5, 1]), 
    np.array([[0,0,0,0],
              [1767732205903/2027836641118,0,0,0],
              [5535828885825/10492691773637,788022342437/10882634858940,0,0],
              [6485989280629/16251701735622, -4246266847089/9704473918619,10755448449292/10357097424841,0]]), 
    np.array([1471266399579/7840856788654, -4482444167858/7529755066697, 11266239266428/11593286722821, 1767732205903/4055673282236]))
tableau_3_i = Tableau(
    np.array([0, 1767732205903/2027836641118, 3/5, 1]), 
    np.array([[0,0,0,0],
              [1767732205903/4055673282236,1767732205903/4055673282236,0,0],
              [2746238789719/10658868560708, -640167445237/6845629431997, 1767732205903/4055673282236, 0],
              [1471266399579/7840856788654, -4482444167858/7529755066697, 11266239266428/11593286722821, 1767732205903/4055673282236]]), 
    np.array([1471266399579/7840856788654, -4482444167858/7529755066697, 11266239266428/11593286722821, 1767732205903/4055673282236]))
tableau_4_e = Tableau(
    np.array([0, 1/2, 83/250, 31/50, 17/20, 1]), 
    np.array([[0,0,0,0,0,0], 
              [1/4, 1/4, 0, 0 ,0 ,0], 
              [8611/62500, -1743/31250 ,1/4 ,0 ,0 ,0], 
              [5012029/34652500, -654441/2922500, 174375/388108, 1/4, 0, 0], 
              [15267082809/155376265600 ,-71443401/120774400 ,730878875/902184768 ,2285395/8070912 ,1/4 ,0], 
              [82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]]), 
    np.array([82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]))
tableau_4_i = Tableau(
    np.array([0, 1/2, 83/250, 31/50, 17/20, 1]), 
    np.array([[0,0,0,0,0,0], 
              [1/4, 1/4, 0, 0 ,0 ,0], 
              [8611/62500, -1743/31250 ,1/4 ,0 ,0 ,0], 
              [5012029/34652500, -654441/2922500, 174375/388108, 1/4, 0, 0], 
              [15267082809/155376265600 ,-71443401/120774400 ,730878875/902184768 ,2285395/8070912 ,1/4 ,0], 
              [82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]]), 
    np.array([82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]))
tableau_5_e = Tableau(
    np.array([0, 0.5, 83/250, 31/50, 17/20, 1]), 
    np.array([[0,0,0,0,0,0],
              [0.5,0,0,0,0,0], 
              [13861/62500, 6889/62500, 0, 0, 0, 0],
              [-116923316275/2393684061468, -2731218467317/15368042101831, 9408046702089/11113171139209, 0, 0, 0],
              [-451086348788/2902428689909, -2682348792572/7519795681897,12662868775082/11960479115383, 3355817975965/11060851509271,0,0],
              [647845179188/3216320057751, 73281519250/8382639484533, 552539513391/3454668386233, 3354512671639/8306763924573, 4040/17871, 0]]), 
    np.array([4586570599.0/29645900160.0, 0, 178811875.0/945068544.0, 814220225.0/1159782912.0, -3700637.0/11593932.0, 61727.0/225920.0]))
tableau_5_i = Tableau(
    np.array([0, 1/2, 83/250, 31/50, 17/20, 1]), 
    np.array([[0,0,0,0,0, 0],
              [1/4,1/4,0,0,0,0],
              [8611/62500, -1743/31250, 1/4, 0, 0, 0],
              [5012029/34652500, -654441/2922500, 174375/388108, 1/4, 0, 0],
              [15267082809/155376265600, -71443401/120774400, 730878875/902184768, 2285395/8070912, 1/4, 0],
              [82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]]), 
    np.array([4586570599.0/29645900160.0, 0, 178811875.0/945068544.0, 814220225.0/1159782912.0, -3700637.0/11593932.0, 61727.0/225920.0]))
tableau_e_embedded = EmbeddedTableau(
    np.array([0, 1/2, 83/250, 31/50, 17/20, 1]), 
    np.array([[0,0,0,0,0,0],[1/2,0,0,0,0,0], 
              [13861/62500, 6889/62500,0,0,0,0], 
              [-116923316275/2393684061468, -2731218467317/15368042101831, 9408046702089/11113171139209, 0, 0, 0], 
              [-451086348788/2902428689909, -2682348792572/7519795681897, 12662868775082/11960479115383, 3355817975965/11060851509271, 0, 0], 
              [647845179188/3216320057751, 73281519250/8382639484533, 552539513391/3454668386233, 3354512671639/8306763924573, 4040/17871 , 0]]), 
    np.array([82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]), 
    np.array([4586570599.0/29645900160.0, 0, 178811875.0/945068544.0, 814220225.0/1159782912.0, -3700637.0/11593932.0, 61727.0/225920.0]), 3)
tableau_i_embedded = EmbeddedTableau(
    np.array([0, 1/2, 83/250, 31/50, 17/20, 1]), 
    np.array([[0,0,0,0,0,0], 
              [1/4, 1/4, 0, 0 ,0 ,0], 
              [8611/62500, -1743/31250 ,1/4 ,0 ,0 ,0], 
              [5012029/34652500, -654441/2922500, 174375/388108, 1/4, 0, 0], 
              [15267082809/155376265600 ,-71443401/120774400 ,730878875/902184768 ,2285395/8070912 ,1/4 ,0], 
              [82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]]), 
    np.array([82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4]), 
    np.array([4586570599.0/29645900160.0, 0, 178811875.0/945068544.0, 814220225.0/1159782912.0, -3700637.0/11593932.0, 61727.0/225920.0]), 3)

class TestExact:
    def test_order_1(self):
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
        assert ark_solve([f1,f2],dt,y0,t0,tf,[tableau_1_e,tableau_1_i]) == yf

    def test_order_3(self):
        A1 = np.array([1,0])
        B1 = -0.5
        C1 = np.array([2,3])
        def f1(t,y):
            return (A1 * t*t + B1*t+C1)
        A2 = np.array([1,5])
        B2 = 1
        C2 = np.array([1,0])
        def f2(t,y):
            return (A2 * t*t + B2*t+C2)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A1+A2)*tf**3/3 + (B1+B2)*tf**2/2 + (C1+C2)*tf + y0)
        assert ark_solve([f1,f2],dt,y0,t0,tf,[tableau_3_e,tableau_3_i]) == yf

    def test_order_4(self):
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
        assert ark_solve([f1,f2],dt,y0,t0,tf,[tableau_4_e,tableau_4_i]) == yf


    def test_order_5(self):
        A1 = np.array([1,0])
        B1 = -0.5
        C1 = np.array([2,3])
        D1 = np.array([-1,2])
        def f1(t,y):
            return (A1 * t**4 + B1*t**3+C1*t+D1)
        A2 = np.array([0,1])
        B2 = 1
        C2 = np.array([1,0])
        D2 = np.array([2,4])
        def f2(t,y):
            return (A2 * t**4 + B2*t**2+C2*t+D2)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A1+A2)*tf**5/5 + (B1)*tf**4/4 + B2*tf**3/3 + (C1+C2)*tf**2/2 + (D1+D2)*tf + y0)
        assert ark_solve([f1,f2],dt,y0,t0,tf,[tableau_5_e,tableau_5_i]) == yf

class TestOrder:
    # Basic Test problem
    def f1(t,y):
        x = y[0]
        y = y[1]
        return np.array([(-y*t)/t, (x*t)/t])
    def f2(t,y):
        x = y[0]
        y = y[1]
        return np.array([(+x)/t, (-y)/t])

    y0 = np.array([np.pi,-np.pi])
    t0 = np.pi
    tf = 2*np.pi
    yf = np.array([-np.pi,1/2+np.pi])

    # non-linear test problem
    B = 20
    C = 4
    def f_non_linear(t,y):
        x = y[0]
        y = y[1]
        return np.array([0,x/TestOrder.B])
    def f2_nl(t,y):
        x = y[0]
        y = y[1]
        return np.array([-x**2/TestOrder.B,0])
    y0_nl = np.array([C,np.log(1/C)])
    t0_nl = 0
    tf_nl = B
    yf_nl = np.array([C/(C+1),np.log((C+1)/C)])


    @pytest.mark.parametrize("method,order,dt0", [([tableau_1_i, tableau_1_e],1,np.pi/8),
                                                  ([tableau_1_e,tableau_1_i],1,np.pi/4),
                                                  ([tableau_3_i,tableau_3_e],3,np.pi/4),
                                                  ([tableau_4_i,tableau_4_e],4,np.pi/4)])
    def test_order(self,method,order,dt0):
        for dt in [dt0,dt0/2,dt0/4]:
            assert np.log2(np.linalg.norm(ark_solve([TestOrder.f1,TestOrder.f2],dt,TestOrder.y0,TestOrder.t0,TestOrder.tf,method)-TestOrder.yf)/(np.linalg.norm(ark_solve([TestOrder.f1,TestOrder.f2],dt/2,TestOrder.y0,TestOrder.t0,TestOrder.tf,method)-TestOrder.yf))) == pytest.approx(order, abs=0.4, rel=0.1)

    @pytest.mark.parametrize("method,order,dt0", [([tableau_1_i, tableau_1_e],1,1/8),
                                                  ([tableau_3_i,tableau_3_e],3,1/4),
                                                  ([tableau_4_i,tableau_4_e],4,tf_nl/2)])
    def test_order_non_linear(self,method,order,dt0):
        for dt in [dt0,dt0/2,dt0/4]:
            assert np.log2(np.linalg.norm(ark_solve([TestOrder.f_non_linear,TestOrder.f2_nl],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,method)-TestOrder.yf_nl)/(np.linalg.norm(ark_solve([TestOrder.f_non_linear,TestOrder.f2_nl],dt/2,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,method)-TestOrder.yf_nl))) == pytest.approx(order, abs=0.4, rel=0.1)
