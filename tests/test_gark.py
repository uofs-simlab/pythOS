import numpy as np
import pytest
import sys
sys.path.insert(1,sys.path[0]+"/../")
from butcher_tableau import Tableau, EmbeddedTableau
from gark_methods import gark_solve

class TestGARK:
    def test_exact(self):
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

        beta = 0.5
        Aee = np.array([[0, 0, 0], [1/2, 0, 0], [1-beta, beta, 0]])
        Aei = np.array([[0, 0], [1/2, 0], [1/2, 1/2]])
        Aie = np.array([[1/4, 0, 0], [1/4, 1/2, 0]])
        Aii = np.array([[1/4, 0], [1/2, 1/4]])

        bE = np.array([1/4, 1/2, 1/4])
        bI = np.array([1/2, 1/2])

        A = [[Aee, Aei], [Aie, Aii]]
        b = [bE, bI]

        assert gark_solve([f1,f2], dt, y0, t0, tf, A, b) == yf

    def test_order(self):
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
        beta = 0.5
        Aee = np.array([[0, 0, 0], [1/2, 0, 0], [1-beta, beta, 0]])
        Aei = np.array([[0, 0], [1/2, 0], [1/2, 1/2]])
        Aie = np.array([[1/4, 0, 0], [1/4, 1/2, 0]])
        Aii = np.array([[1/4, 0], [1/2, 1/4]])

        bE = np.array([1/4, 1/2, 1/4])
        bI = np.array([1/2, 1/2])

        A = [[Aee, Aei], [Aie, Aii]]
        b = [bE, bI]
        dt0 = np.pi/4
        for dt in [dt0,dt0/2,dt0/4]:
            assert np.log2(np.linalg.norm(gark_solve([f1, f2], dt, y0, t0, tf, A, b) - yf) / np.linalg.norm(gark_solve([f1, f2], dt/2, y0, t0, tf, A, b) - yf)) == pytest.approx(2, abs=0.2)
        # non-linear test problem
        B = 1
        C = 2
        def f1(t,y):
            x = y[0]
            y = y[1]
            return np.array([0,x/B])
        def f2(t,y):
            x = y[0]
            y = y[1]
            return np.array([-x**2/B,0])
        y0 = np.array([C,np.log(1/C)])
        t0 = 0
        tf = B
        yf = np.array([C/(C+1),np.log((C+1)/C)])
        dt0 = 1/8
        for dt in [dt0,dt0/2,dt0/4]:
            assert np.log2(np.linalg.norm(gark_solve([f1, f2], dt, y0, t0, tf, A, b) - yf) / np.linalg.norm(gark_solve([f1, f2], dt/2, y0, t0, tf, A, b) - yf)) == pytest.approx(2, abs=0.2)
