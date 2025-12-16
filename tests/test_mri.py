import numpy as np
import math as m
import pytest
import sys
sys.path.insert(1,sys.path[0]+"/../")
from multirate_infinitesimal import *

class TestExact:
    
    @pytest.mark.parametrize("method", [mri_erk2a, mri_erk2b, mri_irk2])
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
        assert (multirate_infinitesimal_solve(y0,t0,dt,tf,method,f1,f2) == yf)

    @pytest.mark.parametrize("method", [mri_kw3, mri_esdirk3a, mri_imex3])
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
        A3 = np.array([2,1])
        B3 = 1
        C3 = np.array([0.5, 3])
        def f3(t, y):
            return (A3 * t*t + B3*t+C3)

        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx((A1+A2)*tf**3/3 + (B1+B2)*tf**2/2 + (C1+C2)*tf + y0)
        yf2 = pytest.approx((A1+A2+A3)*tf**3/3 + (B1+B2+B3)*tf**2/2 + (C1+C2+C3)*tf + y0)

        assert (multirate_infinitesimal_solve(y0,t0,dt,tf,method,f1,f2) == yf)
        if method.omeg is not None:
            assert (multirate_infinitesimal_solve(y0,t0,dt,tf,method,f1,f2, f3) == yf2)

    @pytest.mark.parametrize("method", [mri_erk4, mri_imex4])
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
        yf = pytest.approx((A1+A2)*tf**4/4 + (B1+B2)*tf**3/3 + (C1+C2)*tf**2/2 + (D1+D2)*tf + y0)
        yf2 = pytest.approx((A1+A2+A3)*tf**4/4 + (B1+B2+B3)*tf**3/3 + (C1+C2+C3)*tf**2/2 + (D1+D2+D3)*tf + y0)

        assert (multirate_infinitesimal_solve(y0,t0,dt,tf,method,f1,f2) == yf)
        if method.omeg is not None:
            assert (multirate_infinitesimal_solve(y0,t0,dt,tf,method,f1,f2, f3) == yf2)


class TestOrder:
    lambdaf = -10
    lambdas = -1
    eps = 0.1
    alpha = 1
    beta = 20

    def fe(t, y):
        v = y[1]
        return np.array([0, -np.sin(t)/(2*v)])
    arr =np.array([[lambdaf, (1-eps)/alpha * (lambdaf - lambdas)], [-alpha*eps*(lambdaf-lambdas), lambdas]])

    def fs(t, y):
        u = y[0]
        v = y[1]
        result = TestOrder.arr @ np.array([(-3 + u**2 - np.cos(TestOrder.beta*t))/(2*u), (-2+v**2-np.cos(t))/(2*v)]) - np.array([TestOrder.beta*np.sin(TestOrder.beta*t)/(2*u), np.sin(t)/(2*v)])
        result[0] = 0
        return result

    def ff(t, y):
        u = y[0]
        v = y[1]
        result = TestOrder.arr @ np.array([(-3 + u**2 - np.cos(TestOrder.beta*t))/(2*u), (-2+v**2-np.cos(t))/(2*v)]) - np.array([TestOrder.beta*np.sin(TestOrder.beta*t)/(2*u), np.sin(t)/(2*v)])
        result[1] = 0
        return result

    t0 = 0
    tf = 5*np.pi/2
    y0 = np.array([2, 3**0.5])

    yf = np.array([np.sqrt(3+np.cos(beta*tf)), np.sqrt(2+np.cos(tf))])
    def fi(t, y):
        u = y[0]
        v = y[1]
        result = TestOrder.arr @ np.array([(-3 + u**2 - np.cos(TestOrder.beta*t))/(2*u), (-2+v**2-np.cos(t))/(2*v)])
        result[0] = 0
        return result

    @pytest.mark.parametrize("method,order,dt0,name", [(mri_erk2a, 2, tf/32, "mri_erk2a"),
                                                       (mri_erk2b, 2, tf/32, "mri_erk2b"),
                                                       (mri_irk2, 2, tf/64, "mri_irk2"),
                                                       (mri_kw3, 3, tf/32, "mri_kw3"),
                                                       (mri_esdirk3a, 3, tf/8, "mri_esdirk3a"),
                                                       (mri_imex3, 3, tf/8, "mri_imex3"),
                                                       (mri_erk4, 4, tf/64, "mri_erk4"),
                                                       (mri_imex4, 4, tf/64, "mri_imex4")])
    def test_order(self,method,order,dt0,name):
        for dt in [dt0,dt0/2,dt0/4]:
            assert np.log2(np.linalg.norm(multirate_infinitesimal_solve(TestOrder.y0,TestOrder.t0,dt,TestOrder.tf,method,TestOrder.fs,TestOrder.ff)-TestOrder.yf) / np.linalg.norm(multirate_infinitesimal_solve(TestOrder.y0,TestOrder.t0,dt/2,TestOrder.tf,method,TestOrder.fs,TestOrder.ff)-TestOrder.yf)) == pytest.approx(order, abs=0.4, rel=0.1)

    @pytest.mark.parametrize("method,order,dt0,name", [(mri_imex3, 3, tf/64, "mri_imex3"),
                                                       (mri_imex4, 4, tf/32, "mri_imex4")])
    def test_order_imex(self,method,order,dt0,name):
        for dt in [dt0,dt0/2,dt0/4]:
            assert np.log2(np.linalg.norm(multirate_infinitesimal_solve(TestOrder.y0,TestOrder.t0,dt,TestOrder.tf,method,TestOrder.fi,TestOrder.ff, TestOrder.fe)-TestOrder.yf) / np.linalg.norm(multirate_infinitesimal_solve(TestOrder.y0,TestOrder.t0,dt/2,TestOrder.tf,method,TestOrder.fi,TestOrder.ff, TestOrder.fe)-TestOrder.yf)) == pytest.approx(order, abs=0.4, rel=0.1)
