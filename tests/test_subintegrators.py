import pytest
import numpy as np
import sys
sys.path.insert(1,sys.path[0]+"/../")
import fractional_step as fs

try:
    import sundials
    has_sundials = True
except:
    has_sundials = False

class TestExact:
    @pytest.mark.parametrize("method", ["FE","BE"])
    def test_order_1(self,method):
        A = 2
        def f(t,y):
            return np.ones(y.size)*A
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx(A*tf + y0)
        assert (fs.fractional_step([f],dt,y0,t0,tf,"Godunov",{(0,):method}) == yf)
    
    @pytest.mark.parametrize("method", ["Heun", "GL2", "CN", "SD2O2", "SDAstable","SDLstable","EPI2","EPIRK2"])
    def test_order_2(self,method):
        A = 2
        B = -0.5
        def f(t,y):
            return np.ones(y.size)*(A * t + B)

        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx(A*tf*tf*0.5 + B*tf + y0)
        assert (fs.fractional_step([f],dt,y0,t0,tf,"Godunov",{(0,):method}) == yf)

    @pytest.mark.parametrize("method", ["RK3", "SSPRK3","SD2O3","SD3O3Lstable","EPI3","EPIRK3"])
    def test_order_3(self,method):
        A = 2
        B = -0.5
        C = 1
        def f(t,y):
            return np.ones(y.size)*(A * t*t + B*t+C)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx(A/3*tf**3 + B/2*tf**2 + C*tf + y0)
        assert (fs.fractional_step([f],dt,y0,t0,tf,"Godunov",{(0,):method}) == yf)

    @pytest.mark.parametrize("method", ["RK4","GL4","SSP(5,4)","SD3O4","SD5O4","EPI4","EPIRK4","EPIRK4s3"])
    def test_order_4(self,method):
        A = 2
        B = -0.5
        C = 1
        D = 0.8
        def f(t,y):
            return np.ones(y.size)*(A * t**3 + B*t**2 + C*t + D)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.4
        yf = pytest.approx(A/4*tf**4 + B/3*tf**3 + C/2*tf**2 + D*tf + y0,rel=1e-3)
        assert (fs.fractional_step([f],dt,y0,t0,tf,"Godunov",{(0,):method}) == yf)

    @pytest.mark.parametrize("method", ["EPI5","EPIRK5s3","EPIRK5s4"])
    def test_order_5(self,method):
        A = 2
        B = -0.5
        C = 1
        D = 0.2
        def f(t,y):
            return np.ones(y.size)*(A * t**4 + B*t**3 + C*t + D)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx(A/5*tf**5 + B/4*tf**4 + C/2*tf**2 + D*tf + y0)
        assert (fs.fractional_step([f],dt,y0,t0,tf,"Godunov",{(0,):method}) == yf)

    @pytest.mark.parametrize("method", ["GL6","EPI6"])
    def test_order_6(self,method):
        A = 2
        B = -0.5
        C = 1
        D = 0.2
        def f(t,y):
            return np.ones(y.size)*(A * t**5 + B*t**3 + C*t + D)
        y0 = np.array([0,0.5])
        t0 = 0
        tf = 2
        dt = 0.1
        yf = pytest.approx(A/6*tf**6 + B/4*tf**4 + C/2*tf**2 + D*tf + y0)
        assert (fs.fractional_step([f],dt,y0,t0,tf,"Godunov",{(0,):method}) == yf)

class TestOrder:
    # Basic Test problem
    def f(t,y):
        x = y[0]
        y = y[1]
        return np.array([(-y*t+x)/t, (x*t-y)/t])
    y0 = np.array([np.pi,-np.pi])
    t0 = np.pi
    tf = 2*np.pi
    yf = np.array([-np.pi,1/2+np.pi])

    # non-linear test problem for EPI methods
    B = 20
    C = 4
    def f_non_linear(t,y):
        x = y[0]
        y = y[1]
        return np.array([-x**2/TestOrder.B,x/TestOrder.B])
    y0_nl = np.array([C,np.log(1/C)])
    t0_nl = 0
    tf_nl = B
    yf_nl = np.array([C/(C+1),np.log((C+1)/C)])

    @pytest.mark.parametrize("method,order,dt0", [("FE",1,np.pi/4),
                                                  ("BE",1,np.pi/4),
                                                  ("Heun",2,np.pi/4),
                                                  ("GL2",2,np.pi/4),
                                                  ("CN",2,np.pi/4),
                                                  ("SD2O2",2,np.pi/4),
                                                  ("SDAstable",2,np.pi/4),
                                                  ("SDLstable",2,np.pi/16),
                                                  ("RK3",3,np.pi/4),
                                                  ("SSPRK3",3,np.pi/4),
                                                  ("SD2O3",3,np.pi/4),
                                                  ("SD3O3Lstable",3,np.pi/4),
                                                  ("RK4",4,np.pi/4),
                                                  ("GL4",4,np.pi/4),
                                                  ("SSP(5,4)",4,np.pi/4),
                                                  ("SD3O4",4,np.pi/8),
                                                  ("SD5O4",4,np.pi/4),
                                                  ("GL6",6,np.pi/4)])
    def test_order(self,method,order,dt0):
        for dt in [dt0,dt0/2,dt0/4]:
            assert np.log2(np.linalg.norm(fs.fractional_step([TestOrder.f],dt,TestOrder.y0,TestOrder.t0,TestOrder.tf,"Godunov",{(0,):method})-TestOrder.yf)/(np.linalg.norm(fs.fractional_step([TestOrder.f],dt/2,TestOrder.y0,TestOrder.t0,TestOrder.tf,"Godunov",{(0,):method})-TestOrder.yf))) == pytest.approx(order, abs=0.4, rel=0.1)
    
    # We already have a complete order test of the RK methods, so just run a sampling with the non-linear problem
    @pytest.mark.parametrize("method,order,dt0", [("FE",1,1/4),
                                                  ("BE",1,1/4),
                                                  ("Heun",2,1/4),
                                                  ("GL2",2,1/4),
                                                  ("CN",2,1/4),
                                                  ("SD2O2",2,1/4),
                                                  ("RK3",3,1/4),
                                                  ("SD2O3",3,1),
                                                  ("GL4",4,tf_nl/5),
                                                  ("EPIRK2",2,1/4),
                                                  ("EPIRK3",3,1/4),
                                                  ("EPIRK4s3",4,5),
                                                  ("EPIRK4",4,5),
                                                  ("EPIRK5s3",5,tf_nl/6),
                                                  ("EPIRK5s4",5,tf_nl/5)])
    def test_order_nonlinear(self,method,order,dt0):
        for dt in [dt0,dt0/2,dt0/4]:
            result1 = fs.fractional_step([TestOrder.f_non_linear],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method})
            result2 = fs.fractional_step([TestOrder.f_non_linear],dt/2,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method})
            assert np.log2(np.linalg.norm((result1-TestOrder.yf_nl)/TestOrder.yf_nl)/np.linalg.norm((result2-TestOrder.yf_nl)/TestOrder.yf_nl)) ==pytest.approx(order, abs=0.5, rel=0.1)
    
    @pytest.mark.parametrize("method,order,dt0", [("EPI2",2,tf_nl),
                                                  ("EPI3",3,tf_nl),
                                                  ("EPI4",4,tf_nl)])
    def test_order_nonlinear_epi(self,method,order,dt0):
        for dt in [dt0,dt0/2,dt0/4]:
            result1 = fs.fractional_step([TestOrder.f_non_linear],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters={1:{'monolithic': True}})
            result2 = fs.fractional_step([TestOrder.f_non_linear],dt/2,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters={1:{'monolithic': True}})
            assert np.log2(np.linalg.norm((result1-TestOrder.yf_nl)/TestOrder.yf_nl)/np.linalg.norm((result2-TestOrder.yf_nl)/TestOrder.yf_nl)) ==pytest.approx(order, abs=0.5, rel=0.1)
    
    @pytest.mark.parametrize("method,order,dt0", [("EPI2",2,tf_nl),
                                                  ("EPI3",3,tf_nl),
                                                  ("EPI4",4,tf_nl),
                                                  ("EPIRK2",2,1/4),
                                                  ("EPIRK3",3,1/4),
                                                  ("EPIRK4s3",4,5),
                                                  ("EPIRK4",4,5),
                                                  ("EPIRK5s3",5,tf_nl/6),
                                                  ("EPIRK5s4",5,tf_nl/5)])
    def test_order_nonlinear_epi_exode_scipy(self,method,order,dt0):
        if 'EPIRK' in method:
            solver_parameters={}
        else:
            solver_parameters={1:{'monolithic':True}}
        for dt in [dt0,dt0/2,dt0/4]:
            result1 = fs.fractional_step([TestOrder.f_non_linear],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:('RK45',(1e-10,1e-12))})
            result2 = fs.fractional_step([TestOrder.f_non_linear],dt/2,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:('RK45',(1e-10,1e-12))})
            assert np.log2(np.linalg.norm((result1-TestOrder.yf_nl)/TestOrder.yf_nl)/np.linalg.norm((result2-TestOrder.yf_nl)/TestOrder.yf_nl)) ==pytest.approx(order, abs=0.5, rel=0.1)
    
    @pytest.mark.parametrize("method,order,dt0", [("EPI2",2,tf_nl),
                                                  ("EPI3",3,tf_nl),
                                                  ("EPI4",4,tf_nl),
                                                  ("EPIRK2",2,1/4),
                                                  ("EPIRK3",3,1/4),
                                                  ("EPIRK4s3",4,5),
                                                  ("EPIRK4",4,5),
                                                  ("EPIRK5s3",5,tf_nl/6),
                                                  ("EPIRK5s4",5,tf_nl/5)])
    def test_order_nonlinear_epi_exode_pythOS(self,method,order,dt0):
        if 'EPIRK' in method:
            solver_parameters={}
        else:
            solver_parameters={1:{'monolithic':True}}
        for dt in [dt0,dt0/2,dt0/4]:
            result1 = fs.fractional_step([TestOrder.f_non_linear],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:('Cash-Karp',(1e-10,1e-12))})
            result2 = fs.fractional_step([TestOrder.f_non_linear],dt/2,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:('Cash-Karp',(1e-10,1e-12))})
            assert np.log2(np.linalg.norm((result1-TestOrder.yf_nl)/TestOrder.yf_nl)/np.linalg.norm((result2-TestOrder.yf_nl)/TestOrder.yf_nl)) ==pytest.approx(order, abs=0.5, rel=0.1)

    @pytest.mark.skipif(not has_sundials, reason="requires Sundials dependency")
    @pytest.mark.parametrize("method,order,dt0", [("EPI2",2,tf_nl),
                                                  ("EPI3",3,tf_nl),
                                                  ("EPI4",4,tf_nl),
                                                  ("EPIRK2",2,1/4),
                                                  ("EPIRK3",3,1/4),
                                                  ("EPIRK4s3",4,5),
                                                  ("EPIRK4",4,5),
                                                  ("EPIRK5s3",5,tf_nl/6),
                                                  ("EPIRK5s4",5,tf_nl/5)])
    def test_order_nonlinear_epi_exode_cvode(self,method,order,dt0):
        if 'EPIRK' in method:
            solver_parameters={}
        else:
            solver_parameters={1:{'monolithic':True}}
        for dt in [dt0,dt0/2,dt0/4]:
            result1 = fs.fractional_step([TestOrder.f_non_linear],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:('CV_BDF',(1e-10,1e-12))})
            result2 = fs.fractional_step([TestOrder.f_non_linear],dt/2,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:('CV_BDF',(1e-10,1e-12))})
            assert np.log2(np.linalg.norm((result1-TestOrder.yf_nl)/TestOrder.yf_nl)/np.linalg.norm((result2-TestOrder.yf_nl)/TestOrder.yf_nl)) ==pytest.approx(order, abs=0.5, rel=0.1)

    @pytest.mark.skipif(not has_sundials, reason="requires Sundials dependency")
    @pytest.mark.parametrize("method,order,dt0", [("EPI2",2,tf_nl),
                                                  ("EPI3",3,tf_nl),
                                                  ("EPI4",4,tf_nl),
                                                  ("EPIRK2",2,1/4),
                                                  ("EPIRK3",3,1/4),
                                                  ("EPIRK4s3",4,5),
                                                  ("EPIRK4",4,5),
                                                  ("EPIRK5s3",5,tf_nl/6),
                                                  ("EPIRK5s4",5,tf_nl/5)])
    def test_order_nonlinear_epi_exode_arkode(self,method,order,dt0):
        if 'EPIRK' in method:
            solver_parameters={}
        else:
            solver_parameters={1:{'monolithic':True}}
        for dt in [dt0,dt0/2,dt0/4]:
            result1 = fs.fractional_step([TestOrder.f_non_linear],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:(("ARKODE_ARK548L2SA_DIRK_8_4_5","ARKODE_ERK_NONE"),(1e-10,1e-12))})
            result2 = fs.fractional_step([TestOrder.f_non_linear],dt/2,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:(("ARKODE_ARK548L2SA_DIRK_8_4_5","ARKODE_ERK_NONE"),(1e-10,1e-12))})
            assert np.log2(np.linalg.norm((result1-TestOrder.yf_nl)/TestOrder.yf_nl)/np.linalg.norm((result2-TestOrder.yf_nl)/TestOrder.yf_nl)) ==pytest.approx(order, abs=0.5, rel=0.1)

    @pytest.mark.skipif(not has_sundials, reason="requires Sundials dependency")
    @pytest.mark.parametrize("method,order,dt0", [("EPI2",2,tf_nl),
                                                  ("EPI3",3,tf_nl),
                                                  ("EPI4",4,tf_nl),
                                                  ("EPIRK2",2,1/4),
                                                  ("EPIRK3",3,1/4),
                                                  ("EPIRK4s3",4,5),
                                                  ("EPIRK4",4,5),
                                                  ("EPIRK5s3",5,tf_nl/6),
                                                  ("EPIRK5s4",5,tf_nl/5)])
    def test_order_nonlinear_epi_exode_erk_step(self,method,order,dt0):
        if 'EPIRK' in method:
            solver_parameters={}
        else:
            solver_parameters={1:{'monolithic':True}}
        for dt in [dt0,dt0/2,dt0/4]:
            result1 = fs.fractional_step([TestOrder.f_non_linear],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:('ARKODE_DORMAND_PRINCE_7_4_5',(1e-10,1e-12))})
            result2 = fs.fractional_step([TestOrder.f_non_linear],dt/2,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):method},solver_parameters=solver_parameters,epi_options={1:('ARKODE_DORMAND_PRINCE_7_4_5',(1e-10,1e-12))})
            assert np.log2(np.linalg.norm((result1-TestOrder.yf_nl)/TestOrder.yf_nl)/np.linalg.norm((result2-TestOrder.yf_nl)/TestOrder.yf_nl)) ==pytest.approx(order, abs=0.5, rel=0.1)

class TestEmbeddedTableau:
    # Basic Test problem
    def f(t,y):
        x = y[0]
        y = y[1]
        return np.array([(-y*t+x)/t, (x*t-y)/t])
    y0 = np.array([np.pi,-np.pi])
    t0 = np.pi
    tf = 2*np.pi
    yf = np.array([-np.pi,1/2+np.pi])

    # non-linear test problem for EPI methods
    B = 20
    C = 4
    def f_non_linear(t,y):
        x = y[0]
        y = y[1]
        return np.array([-x**2/TestOrder.B,x/TestOrder.B])
    y0_nl = np.array([C,np.log(1/C)])
    t0_nl = 0
    tf_nl = B
    yf_nl = np.array([C/(C+1),np.log((C+1)/C)])
    @pytest.mark.parametrize("method", ["Dormand-Prince","Cash-Karp","Fehlberg","Bogacki-Shampine","Heun-Euler"])
    def test_embedded(self,method):
        for dt in [TestEmbeddedTableau.tf_nl]:
            result1 = fs.fractional_step([TestOrder.f],dt,TestOrder.y0,TestOrder.t0,TestOrder.tf,"Godunov",{(0,):"ADAPTIVE"},ivp_methods={1:(method,1e-10,1e-12)})
            assert result1 == pytest.approx(TestOrder.yf,rel=3e-5)

    @pytest.mark.parametrize("method", ["Dormand-Prince","Cash-Karp","Fehlberg","Bogacki-Shampine","Heun-Euler"])
    def test_nonlinear_embedded(self,method):
        for dt in [TestEmbeddedTableau.tf_nl]:
            result1 = fs.fractional_step([TestOrder.f_non_linear],dt,TestOrder.y0_nl,TestOrder.t0_nl,TestOrder.tf_nl,"Godunov",{(0,):"ADAPTIVE"},ivp_methods={1:(method,1e-10,1e-12)})
            assert result1 == pytest.approx(TestOrder.yf_nl,rel=1e-5)
