#include "gas.h"
#include <cmath>
#include <iostream>
#include <fstream>

// constructor
gas::gas(controlProperties *contProp_, systemProperties *sysProp_, bool condAW, bool trumpetlobule, int sp, int wo, double gridsize, double length){

    contProp = contProp_;
    sysProp  = sysProp_;

    double pi = 4.*atan(1.);

    species = sp;
    washout = wo;

    // spatial discretisation
    if (gridsize > length){
        cout << "ERROR: grid space is bigger than lenght of duct" << endl;
    }

    if (trumpetlobule){
        N = contProp->NxLob;
    }
    else{
        N = int(length/gridsize);
    }

    if (N < contProp->NxDuctMin){
        //cout << "DANGER: less than "<< contProp->NxDuctMin <<" nodes in this duct! " << N << " -> kept at five nodes." << endl;
        N = contProp->NxDuctMin;
    }
    if (N > contProp->NxDuctMax){
        N = contProp->NxDuctMax;
    }

    h = length/double(N);

    l = length;

    diam = 0;

    // Crank-Nicolson coefficient for time integration: theta = 1 => Euler implicite
    theta = contProp->generalCN;

    /* INITIALIZE
     -------------------*/
    u = Q = Deff = Diff = Pe = 0.;

    S = St = 0;

    F_l = F_r = 0;

    kappa      = sysProp->kappa;
    kappa_hat = 2.*kappa*kappa;

    Diff_fac   = sysProp->Diff_fac;

    pParent = pMinorDaughter = pMajorDaughter = pSister = 0;

    // initialize concentration variables
    if (wo == 0){
        if (condAW){
            C     = VectorXd::Ones(N+1);
            cCAW = 1.;
            dcdxCAW = 0.;
        }
        if (trumpetlobule){
            C     = VectorXd::Ones(N+1);
            c1A = c2A = 1.;
        }
    }
    else {
        if (condAW){
            C     = VectorXd::Zero(N+1);
            cCAW = 0.;
            dcdxCAW = 0.;
        }
        if (trumpetlobule){
            C     = VectorXd::Zero(N+1);
            c1A = c2A = 0.;
        }
    }

    // linear system for trumpet power law coefficients
    powLawRhs = VectorXd::Zero(2);
    powLawP   = VectorXd::Zero(2);
    p1 = p2 = p1d = p2d = 0.;
    n1 = sysProp->n1;
    n2 = sysProp->n2;
    z_star = sysProp->z_star;

    powLawM   = MatrixXd::Zero(2,2);

    // for quadratic interpolation
    quadIntA = VectorXd::Zero(3);
    quadIntC = VectorXd::Zero(3);
    quadIntX = VectorXd::Zero(3);

    quadIntM = MatrixXd::Zero(3,3);

    // preconditionning of iterative solver
    solverBiCGSTAB.preconditioner().setDroptol(.001);



    /* WRITE COEFFICIENTS IN FINITE DIFFERENCE OPERATORS (only for conductive airways)
     -------------------*/

    if (condAW){
        chebychevcollocation = false;

        // identity matrices
        I        = MatrixXd::Zero(N+1,N+1);

        // vector with old concentration values
        Cm1 = C;

        // initialize finite difference schemes
        D1_FW    = MatrixXd::Zero(N+1,N+1);
        D1_BW    = MatrixXd::Zero(N+1,N+1);
        D1_C     = MatrixXd::Zero(N+1,N+1);
        D2_C     = MatrixXd::Zero(N+1,N+1);

        // schemes for chebychev collocation
        if (chebychevcollocation){
            D1_CC      = MatrixXd::Zero(N+1,N+1);
            D2_CC      = MatrixXd::Zero(N+1,N+1);
        }

        for (int i=0; i<=N; i++){
            for (int j=0; j<=N; j++){
                // unity matrix
                if (i==j){
                    I(i,j)       = 1.;
                }

                // --- 2nd order ---

                // schemes for first derivatives:
                if (i==j && j>1){
                    D1_BW(i,j)   =  3./2.;
                    D1_BW(i,j-1) = -4./2.;
                    D1_BW(i,j-2) =  1./2.;
                }
                if (i==j && j==1){
                    D1_BW(i,j)   =  3./2.;
                    D1_BW(i,j-1) = -4./2.;
                }
                if (i==j && j==0){
                    D1_BW(i,j)   =  3./2;
                }
                if (i==j && j<N-1){
                    D1_FW(i,j)   = -3./2.;
                    D1_FW(i,j+1) =  4./2.;
                    D1_FW(i,j+2) = -1./2.;
                }
                if (i==j && j==N-1){
                    D1_FW(i,j)   = -3./2.;
                    D1_FW(i,j+1) =  4./2.;
                }
                if (i==j && j==N){
                    D1_FW(i,j)   = -3./2.;
                }

                // schemes for second derivatives
                if (i==j && j>0 && j<N){
                    D2_C(i,j)        = -2.;
                    D2_C(i,j-1)      =  1.;
                    D2_C(i,j+1)      =  1.;
                }
                if (i==j && i==0){
                    D2_C(i,j)        = -2.;
                    D2_C(i,j+1)      =  1.;
                }
                if (i==j && i==N){
                    D2_C(i,j)        = -2.;
                    D2_C(i,j-1)      =  1.;
                }


                // Chenychev-Collocation Derivative Operator
                if (chebychevcollocation){
                    x_i = cos(i*pi/N);
                    x_j = cos(j*pi/N);

                    if (i == 0 || i == N){
                        zeta_i = 2.;
                    }
                    else{
                        zeta_i = 1.;
                    }
                    if (j == 0 || j == N){
                        zeta_j = 2.;
                    }
                    else{
                        zeta_j = 1.;
                    }

                    if (i!=j){
                        D1_CC(i,j) = zeta_i*pow((-1.),(i+j))/(zeta_j*(x_i - x_j));
                    }

                    if (i==j && (i!=0 && i!=N)){
                        D1_CC(i,j) = -x_i/(2.*(1 - (x_i*x_i)));
                    }

                    if (i==j && i==0){
                        D1_CC(i,j) = (2.*N*N + 1.)/(6.);
                    }

                    if (i==j && i==N){
                        D1_CC(i,j) = -(2.*N*N + 1.)/(6.);
                    }
                }
            }
        }

        if (chebychevcollocation){
            D1_CC *= 2./l;
            D2_CC = D1_CC*D1_CC;

            // flip up-down and left-right ('reverse' in EIGEN) for more confortable handling of Chebychev-Collocation operators
            D1_CC.reverseInPlace();
            D2_CC.reverseInPlace();

            // delete first and last row. the corresponding entries in the solution vector are given from the node cell and the b.c.
            D1_CC.row(0).setZero();
            D2_CC.row(N).setZero();
        }

    }

    if (trumpetlobule){

        // transport cross-section initialization
        STrd     = VectorXd::Zero(N+1);

        // cross-section initialization
        ATr      = VectorXd::Zero(N+1);

        // per-generation diameter initialization
        dTr      = VectorXd::Zero(N+1);

        // physical velocity initialization
        UphTr    = VectorXd::Zero(N+1);

        // advection and diffusion operators
        UTr      = MatrixXd::Zero(N+1, N+1);
        DIFFTr   = MatrixXd::Zero(N+1, N+1);

        // identity
        I        = MatrixXd::Zero(N+1,N+1);

        // finite difference schemes initialization
        D1_C     = MatrixXd::Zero(N+1,N+1);
        D2_C     = MatrixXd::Zero(N+1,N+1);

        for (int i=0; i<=N; i++){
            for (int j=0; j<=N; j++){

                // identity
                if(i==j)
                    I(i,j)      =  1.;

                // for first derivatives
                if (i==j && (j>0 && j<N)){

                    D1_C(i,j-1) = -1./2.;
                    D1_C(i,j+1) =  1./2.;

                }

                // for second derivatives
                if (i==j && (j>0 && j<N)){

                    D2_C(i,j-1) =  1./1.;
                    D2_C(i,j)   = -2./1.;
                    D2_C(i,j+1) =  1./1.;

                }
            }
        }
    }

    // linear system for quadratic interpolation in ghost cells
    for (int i=0; i<3; i++){

        // matrix defined by controlpoints for quadratic interpolation
        quadIntM(i,0) = (i*h)*(i*h);
        quadIntM(i,1) = (i*h);
        quadIntM(i,2) = 1.;

        quadIntX(i) = (i*h);
    }

    /* READ DIFFUSION COEFFICIENT
     -------------------*/
    // carrier species is O_2
    double csp = 0;
    if (csp == species){
        cout << "ERROR: carrier gas and tracer gas are the same" << endl;
    }

    // number of specie in sheet
    int Nsp = 5;

    // read from sheet
    ifstream inFile;
    DmolMat = MatrixXd::Zero(Nsp,Nsp);
    inFile.open("data/diffCoeff");

    if (inFile.is_open()){
        for (int i=0; i<Nsp; i++){
            for (int j=0; j<Nsp; j++){
                inFile >> DmolMat(i,j);
            }
        }
        inFile.close();
    }

    Dmol = DmolMat(csp,species)*1e-04;
    Dmol *= Diff_fac;
}


gas::~gas(){

}


void gas::setCrossSection(double d){

    double pi = 4.*atan(1.);

    if (d <= 0.){
        cout << "ERROR: airway has no diameter assigned" << endl;
    }

    S = St = (d*d)/4.*pi;
    diam = d;

}


void gas::setTrumpetLobuleTransportDomain(double terminalDuctLength, double VLb0){

    double pi = 4.*atan(1.);

    // declarations
    double x, z, nAW;

    // set length of last conductive airway
    lt = terminalDuctLength;

    // set initial volume of dynamic and 'advective' trumpet lobule
    // the size and geometry of the 'advective' is defined using a fraction of the initial trumpet volume
    VTr = VLb0;
    VTrd = 0.5*VLb0;

    // define cross-section interaction of power law and exponential law (constraint for power law)
    x_star = 0.;
    for (int k=0; k<z_star; k++){
        x_star += lt*pow(kappa,k+1);
    }

    S_star = St*pow(kappa_hat,z_star);

    // solve linear system for power law coefficients
    powLawM(0,0) = pow(x_star, n1);
    powLawM(0,1) = pow(x_star, n2);
    powLawM(1,0) = pow(l, n1+1)/(n1+1);
    powLawM(1,1) = pow(l, n2+1)/(n2+1);

    powLawRhs(0) = S_star - St;

    // ...for dynamic trumpet lobule...
    powLawRhs(1) = VTr - St*l;
    powLawP = powLawM.householderQr().solve(powLawRhs);
    p1 = powLawP(0);
    p2 = powLawP(1);

    // ..., and for 'advective' trumpet lobule
    powLawRhs(1) = VTrd - St*l;
    powLawP = powLawM.householderQr().solve(powLawRhs);
    p1d = powLawP(0);
    p2d = powLawP(1);


    // check whether coefficients are negative
    if ((p1<0) || (p2<0) || (p1d<0) || (p2d<0)){
        cout << "ERROR: Coefficients in power law are negative. Check constraints." << endl;
        cout << "p1 = " << p1 << "; p2 = " << p2 << "; p1d = " << p1d << "; p2d = " << p2d << endl;
        cout << "n1 = " << n1 << "; n2 = " << n2 << "; VTr = " << VTr << "; lTr = " << l << endl;
        cout << "x_star = " << x_star << "; S_star = " << S_star << "St = " << St << "; kappa_hat = " << kappa_hat << endl;
    }

    // set 'advective' cross-section
    for (int j=0; j<=N; j++){

        // streamwise coordinate
        x = j*h;

        // advective cross-section from power law
        STrd(j) = p1d*pow(x,n1) + p2d*pow(x, n2) + St;

        // per generation diameter
        z = log(x*(kappa-1.)/(kappa*lt) + 1.)/log(kappa);
        //d_z = sqrt(4.*St/pi)*pow(kappa,z);
        nAW = pow(2, z);
        dTr(j) = sqrt(4.*STrd(j)/(pi*nAW)); // this is not really a diameter but a quantity for the computation of Diff_z

    }
}


void gas::updateTrumpetLobuleProperties(double VLb){

    double pi = 4.*atan(1.);

    // volume
    VTr = VLb;

    // time depending power law coefficient
    p1 = (n1+1)/pow(l,n1+1) * (VTr - pow(l, n2+1)/(n2+1)*p2 - St*l);

}


void gas::computeEffectiveDiffCoeff(double d){

    if (d <= 0.){
        cout << "ERROR: airway has no lenght assigned" << endl;
        return;
    }

    if (Dmol <= 0.){
        cout << "ERROR: molecular diffusion coefficient not given" << endl;
        return;
    }

    // compute Peclet number
    Pe = d*abs(u)/Dmol;

    // effective Diffusion based on Taylor dispersion theory
    Deff = min(Dmol*(1. + 1./192.*pow(Pe,2)), 1000*Dmol);
    //Deff = Dmol*(1. + 1./192.*pow(Pe,2));

}


void gas::setDiffusionCoefficient(int opt){

    if (Dmol <= 0){
        cout << "ERROR: molecular diffusion coefficient not given" << endl;
        return;
    }

    if (opt == 0){
        Diff = Dmol;
        return;
    }
    if (opt == 1){
        Diff = Deff;
        return;
    }
}


void gas::updateConcentrationInTrumpetLobule(double dt,  double cLastCondAW, double dcdxLastCondAW, double hLastCondAW){

    double pi = 4.*atan(1.);

    double Erel, Eabs;

    if (h <= 0){
        cout << "ERROR: grid spacing is not given" << endl;
    }

    // declarations
    double dSdx, dDiffdx, STr, QTr, uTr, nAW;
    double Diff_p1, Diff_m1, QTr_p1, QTr_m1;
    double Pe_x, Diff_x;

    // write advection and diffusion matrices
    for (int j=0; j<=N; j++){
        for(int i=0; i<=N; i++){

            if (i==j){

                // streamwise coordinate
                x = j*h;

                // derivative of transport cross-section
                dSdx = n1*p1*pow(x, n1-1) + n2*p2*pow(x, n2-1);

                // trumpet cross-section
                STr = p1*pow(x, n1) + p2*pow(x, n2) + St;
                ATr(j) = STr;

                // flow and velocity
                QTr = Q*(1. - pow(x, n1+1)/pow(l, n1+1));
                uTr = QTr/STrd(j);

                // physical velocity u = u(x,t)
                UphTr(j) = uTr;

                // calculate effective diffusion coefficient (Taylor dispersion)
                Pe_x = uTr*dTr(j)/Dmol;
                Diff_x = 1.0 * Dmol*(1 + 1./192.*Pe_x*Pe_x);

                // finite difference estimation of diffusivity gradient
                if (j==0){
                    QTr_p1 = Q*(1. - pow(x+h, n1+1)/pow(l, n1+1));
                    Diff_p1 = 1.0 * Dmol*(1 + 1./192.*pow(QTr_p1/STrd(j+1)*dTr(j+1)/Dmol,2));
                    dDiffdx = (Diff_p1 - Diff_x)/h;
                }
                else if (j==N){
                    QTr_m1 = Q*(1. - pow(x-h, n1+1)/pow(l, n1+1));
                    Diff_m1 = 1.0 * Dmol*(1 + 1./192.*pow(QTr_m1/STrd(j-1)*dTr(j-1)/Dmol,2));
                    dDiffdx = (Diff_x - Diff_m1)/h;
                }
                else{
                    QTr_p1 = Q*(1. - pow(x+h, n1+1)/pow(l, n1+1));
                    QTr_m1 = Q*(1. - pow(x-h, n1+1)/pow(l, n1+1));
                    Diff_p1 = 1.0 * Dmol*(1 + 1./192.*pow(QTr_p1/STrd(j+1)*dTr(j+1)/Dmol,2));
                    Diff_m1 = 1.0 * Dmol*(1 + 1./192.*pow(QTr_m1/STrd(j-1)*dTr(j-1)/Dmol,2));
                    dDiffdx = (Diff_p1 - Diff_m1)/(2*h);
                }

                // effective advection and effective diffusion
                //dDiffdx = 0.;
                UTr(i,j)    = (STrd(j)/STr*uTr - Diff_x/STr*dSdx - dDiffdx);
                DIFFTr(i,j) = Diff_x;

                if ((j==2) && promptFlag){
                    //cout << "z = " << z << "; d_z/d = " << d_z/diam << "; Pe = " << Pe << "; Pe_z = " << Pe_z << "; D_eff/D_mol = " << Deff/Dmol << "; D_eff_z/D_mol = " << Diff_z/Dmol << endl;
                }
            }
        }
    }

    // left and right operators
    LO = I - dt*    theta *(-1./h*UTr*D1_C + 1./(h*h)*DIFFTr*D2_C);
    RO = I + dt*(1.-theta)*(-1./h*UTr*D1_C + 1./(h*h)*DIFFTr*D2_C);

    // compute right-hand side
    rhs = RO*C;

    // implement boundary conditions (no directional dependency given)

    // continuity between last conductive airway and trumpet lobule
    LO.row(0).setZero();
    LO(0,0) = 1.;
    rhs(0) = cLastCondAW;

    // smoothness between last conductive airway and trumpet lobule
    /*
    LO.row(1).setZero();
    LO(1,0) = -3./(2.*h); LO(1,1) = 4./(2.*h); LO(1,2) = -1./(2.*h);
    rhs(1) = dcdxLastCondAW;
    */

    // zero gradient at end of trumpet lobule
    LO.row(N).setZero();
    LO(N,N) = 3.; LO(N,N-1) = -4.; LO(N,N-2) = 1.;
    rhs(N) = 0.;

    // solve for concentration:
    //C = LO.householderQr().solve(rhs);

    C = solverBiCGSTAB.compute((LO).sparseView()).solve(rhs);

    // mass-flux at the end of trumpet
    F_l = UTr(0,0)*STrd(0)*C(0) - STrd(0)*DIFFTr(0,0)/(2.*h)*(-3.*C(0) + 4.*C(1)   - 1.*C(2));
    F_r = UTr(N,N)*STrd(N)*C(N) - STrd(N)*DIFFTr(N,N)/(2.*h)*( 3.*C(N) - 4.*C(N-1) + 1.*C(N-2));

    // ghost-cell concentration values computed with quadratic interpolation
    c1A = quadraticInterpolation(1.*hLastCondAW, 0, 0);
    c2A = quadraticInterpolation(2.*hLastCondAW, 0, 0);

    if (promptFlag){
        //cout << "c1A/C(1) = " << c1A/C(1) << "; c2A/C(2) = " << c2A/C(2) << "; c1A, c2A = " << c1A << ", " << c2A << "; hCAW = " << hLastCondAW*1000. << "; x0, x1, x2 = " << quadIntX(0)*1000. << ", " << quadIntX(1)*1000. << ", " << quadIntX(2)*1000. << endl;
    }

    // error check
    Erel = (LO*C - rhs).norm() / rhs.norm();
    Eabs = (LO*C - rhs).norm();

    if (Erel > 0.1){
        cout << "OOPS: relative error in concentration update is bigger than 10 %" << endl;
    }

    // clean up
    UTr.setZero();
    DIFFTr.setZero();
}


void gas::updateConcentrationInDuct(bool inletTrachea, bool reachedEnd, double dt,  double cin, double c1Lobule, double c2Lobule){

    if (h <= 0){
        cout << "ERROR: grid spacing is not given" << endl;
    }

    // declarations
    bool diffDom;
    double Erel, Eabs;
    double c_Agcl1, c_Agcl2, c_Agcr1, c_Agcr2;
    double c_Dgcl,  c_Dgcr;
    double wQs, wQp, wQmaj, wQmin;
    double wSs, wSp, wSmaj, wSmin;
    double c_min01, c_min02, c_maj01, c_maj02, c_p1N, c_p2N, c_s01, c_s02;

    // pick correct upwind scheme
    if (u >= 0){
        if (inletTrachea){
            D1_BW.block(0,0,2,N+1) *= 0.;
            D1_BW(1,0) = -1.;
            D1_BW(1,1) =  1.;
        }
        D1_UW = D1_BW;
    }

    if (u < 0){
        D1_UW = D1_FW;
    }

    // check Peclet number (used later for diffusion-ghost cell definition)
    if (Pe < 1.0){
        // diffusion dominant
        diffDom = true;
    }
    else{
        // advection dominant
        diffDom = false;
    }

    // left and right operators
    LO = I - dt*    theta *(-u/h*D1_UW + Diff/(h*h)*D2_C);
    RO = I + dt*(1.-theta)*(-u/h*D1_UW + Diff/(h*h)*D2_C);

    // compute right-hand side
    rhs = RO*C;

    // BOUNDARY CONDITIONS - for trachea
    if (inletTrachea){
        // ------------------ define ghost cell values -------------------
        // define weights for advection
        wQmaj = abs(pMajorDaughter->Q)/(abs(pMajorDaughter->Q) + abs(pMinorDaughter->Q));
        wQmin = abs(pMinorDaughter->Q)/(abs(pMajorDaughter->Q) + abs(pMinorDaughter->Q));

        // define weights for diffusion
        wSmaj = pMajorDaughter->S/(pMajorDaughter->S + pMinorDaughter->S);
        wSmin = pMinorDaughter->S/(pMajorDaughter->S + pMinorDaughter->S);

        // interpolate concentration in neighoring domains
        c_maj01 = pMajorDaughter->quadraticInterpolation(1.*h, 0, 1); c_min01 = pMinorDaughter->quadraticInterpolation(1.*h, 0, 1);
        c_maj02 = pMajorDaughter->quadraticInterpolation(2.*h, 0, 1); c_min02 = pMinorDaughter->quadraticInterpolation(2.*h, 0, 1);

        // define ghost cell values for advection
        if ((pMajorDaughter->u) < 0 && (pMinorDaughter->u < 0)){ // two feeder
            c_Agcr1 = wQmaj*c_maj01 + wQmin*c_min01;
            c_Agcr2 = wQmaj*c_maj02 + wQmin*c_min02;
        }
        if ((pMajorDaughter->u) < 0 && (pMinorDaughter->u >= 0)){ // one feeder - major daughter
            c_Agcr1 = c_maj01;
            c_Agcr2 = c_maj02;
        }
        if ((pMajorDaughter->u) >= 0 && (pMinorDaughter->u < 0)){ // one feeder - minor daughter
            c_Agcr1 = c_min01;
            c_Agcr2 = c_min02;
        }

        if ((u>0) && (pMajorDaughter->u<0) && (pMinorDaughter->u<0)){ // sink apparent
            c_Agcr1 = c_Agcr2 = Cm1(N);
            //             cout << "NON-PHYSICAL: sink apparent at node" << endl;
        }
        if ((u<0) && (pMajorDaughter->u>0) && (pMinorDaughter->u>0)){ // source apparent
            c_Agcr1 = c_Agcr2 = Cm1(N);
            //             cout << "NON-PHYSICAL: source apparent at node" << endl;
        }

        // define ghost cell values for diffusion
        /*
        if (diffDom){
            c_Dgcr = wSmaj*c_maj01 + wSmin*c_min01;
        }
        else{
            if (u >= 0){
                c_Dgcr = wQmaj*c_maj01 + wQmin*c_min01;
            }
            else{
                c_Dgcr = c_Agcr1;
            }
        }
        */
        c_Dgcr = wSmaj*c_maj01 + wSmin*c_min01;

        // ---------------------------------------------------------------

        if (u >= 0){

            // left b.c.
            LO.block(0,0,1,N+1) *= 0.;
            LO(0,0) = 1.;
            rhs(0)  = cin;

            rhs(1) -= dt*u/(2.*h)*(1.*cin);

            // right b.c.
            rhs(N) += dt*Diff/(h*h)*(1.*c_Dgcr);
        }

        if (u < 0){

            // left b.c.
            LO.block(0,0,1,N+1) *= 0.;
            LO(0,0) = -3.; LO(0,1) = 4.; LO(0,2) = -1.;
            rhs(0) *= 0.;

            // right b.c.
            rhs(N) += dt*Diff/(h*h)*(1.*c_Dgcr);
            rhs(N) -= dt*u/(2.*h)*(4.*c_Agcr1 - 1.*c_Agcr2);

            rhs(N-1) -= dt*u/(2.*h)*(-1.*c_Agcr1);
        }
    }

    // BOUNDARY CONDITIONS - for ending airways
    else if (reachedEnd){
        // ------------------ define ghost cell values -------------------
        // define weights for advection
        wQp   = abs(pParent->Q)/(abs(pSister->Q) + abs(pParent->Q));
        wQs   = abs(pSister->Q)/(abs(pSister->Q) + abs(pParent->Q));

        // define weights for diffusion
        wSp   = pParent->S/(pSister->S + pParent->S);
        wSs   = pSister->S/(pSister->S + pParent->S);

        // interpolate concentration in neighoring domains
        c_p1N   = pParent->quadraticInterpolation(1.*h, 1, 1);  c_s01   = pSister->quadraticInterpolation(1.*h, 0, 1);
        c_p2N   = pParent->quadraticInterpolation(2.*h, 1, 1);  c_s02   = pSister->quadraticInterpolation(2.*h, 0, 1);

        // define ghost cell values for advection
        if ((pParent->u) >= 0 && (pSister->u < 0)){ // two feeder
            c_Agcl1 = wQp*c_p1N + wQs*c_s01;
            c_Agcl2 = wQp*c_p2N + wQs*c_s02;
        }
        if ((pParent->u) >= 0 && (pSister->u >= 0)){ // one feeder - parent
            c_Agcl1 = c_p1N;
            c_Agcl2 = c_p2N;
        }
        if ((pParent->u) < 0 && (pSister->u < 0)){ // one feeder - sister
            c_Agcl1 = c_s01;
            c_Agcl2 = c_s02;
        }

        if ((u<0) && (pSister->u<0) && (pParent->u>0)){ // sink apparent
            c_Agcl1 = c_Agcl2 = Cm1(0);
            //             cout << "NON-PHYSICAL: sink apparent at node" << endl;
        }
        if ((u>0) && (pSister->u>0) && (pParent->u<0)){ // source apparent
            c_Agcl1 = c_Agcl2 = Cm1(0);
            //             cout << "NON-PHYSICAL: source apparent at node" << endl;
        }

        // define ghost cell values for diffusion
        /*
        if (diffDom){
            c_Dgcl = wSp*c_p1N     + wSs*c_s01;
        }
        else{
            if (u >= 0){
                c_Dgcl = c_Agcl1;
            }
            else{
                c_Dgcl = wQp*c_p1N     + wQs*c_s01;
            }
        }
        */

        c_Dgcl = wSp*c_p1N     + wSs*c_s01;

        // ---------------------------------------------------------------


        if (u >= 0){

            // left b.c.
            rhs(0) += dt*Diff/(h*h)*(1.*c_Dgcl);
            rhs(0) -= dt*u/(2.*h)*(-4.*c_Agcl1 + 1.*c_Agcl2);

            rhs(1) -= dt*u/(2.*h)*(1.*c_Agcl1);

            // right b.c.
            rhs(N) += dt*Diff/(h*h)*(1.*c1Lobule);
        }

        if (u < 0){

            // left b.c.
            rhs(0) += dt*Diff/(h*h)*(1.*c_Dgcl);

            // right b.c.
            rhs(N) += dt*Diff/(h*h)*(1.*c1Lobule);
            rhs(N) -= dt*u/(2.*h)*(4.*c1Lobule - 1.*c2Lobule);

            rhs(N-1) -= dt*u/(2.*h)*(-1.*c1Lobule);
        }
    }

    // BOUNDARY CONDITIONS - for airways inside domain
    else{
        // ------------------ define ghost cell values -------------------
        // define weights for advection
        wQp   = abs(pParent->Q)/(abs(pSister->Q) + abs(pParent->Q));
        wQs   = abs(pSister->Q)/(abs(pSister->Q) + abs(pParent->Q));
        wQmaj = abs(pMajorDaughter->Q)/(abs(pMajorDaughter->Q) + abs(pMinorDaughter->Q));
        wQmin = abs(pMinorDaughter->Q)/(abs(pMajorDaughter->Q) + abs(pMinorDaughter->Q));

        // define weights for diffusion
        wSp   = pParent->S/(pSister->S + pParent->S);
        wSs   = pSister->S/(pSister->S + pParent->S);
        wSmaj = pMajorDaughter->S/(pMajorDaughter->S + pMinorDaughter->S);
        wSmin = pMinorDaughter->S/(pMajorDaughter->S + pMinorDaughter->S);

        // interpolate concentration in neighoring domains
        c_p1N   = pParent->quadraticInterpolation(1.*h, 1, 1);  c_s01   = pSister->quadraticInterpolation(1.*h, 0, 1);
        c_p2N   = pParent->quadraticInterpolation(2.*h, 1, 1);  c_s02   = pSister->quadraticInterpolation(2.*h, 0, 1);

        c_maj01 = pMajorDaughter->quadraticInterpolation(1.*h, 0, 1); c_min01 = pMinorDaughter->quadraticInterpolation(1.*h, 0, 1);
        c_maj02 = pMajorDaughter->quadraticInterpolation(2.*h, 0, 1); c_min02 = pMinorDaughter->quadraticInterpolation(2.*h, 0, 1);

        // define ghost cell values for advection
        if ((pParent->u) >= 0 && (pSister->u < 0)){ // two feeder
            c_Agcl1 = wQp*c_p1N + wQs*c_s01;
            c_Agcl2 = wQp*c_p2N + wQs*c_s02;
        }
        if ((pParent->u) >= 0 && (pSister->u >= 0)){ // one feeder - parent
            c_Agcl1 = c_p1N;
            c_Agcl2 = c_p2N;
        }
        if ((pParent->u) < 0 && (pSister->u < 0)){ // one feeder - sister
            c_Agcl1 = c_s01;
            c_Agcl2 = c_s02;
        }

        if ((pMajorDaughter->u) < 0 && (pMinorDaughter->u < 0)){ // two feeder
            c_Agcr1 = wQmaj*c_maj01 + wQmin*c_min01;
            c_Agcr2 = wQmaj*c_maj02 + wQmin*c_min02;
        }
        if ((pMajorDaughter->u) < 0 && (pMinorDaughter->u >= 0)){ // one feeder - major daughter
            c_Agcr1 = c_maj01;
            c_Agcr2 = c_maj02;
        }
        if ((pMajorDaughter->u) >= 0 && (pMinorDaughter->u < 0)){ // one feeder - minor daughter
            c_Agcr1 = c_min01;
            c_Agcr2 = c_min02;
        }

        if ((u>0) && (pMajorDaughter->u<0) && (pMinorDaughter->u<0)){ // sink apparent
            c_Agcr1 = c_Agcr2 = Cm1(N);
            //             cout << "NON-PHYSICAL: sink apparent at node" << endl;
        }
        if ((u<0) && (pMajorDaughter->u>0) && (pMinorDaughter->u>0)){ // source apparent
            c_Agcr1 = c_Agcr2 = Cm1(N);
            //             cout << "NON-PHYSICAL: source apparent at node" << endl;
        }
        if ((u<0) && (pSister->u<0) && (pParent->u>0)){ // sink apparent
            c_Agcl1 = c_Agcl2 = Cm1(0);
            //             cout << "NON-PHYSICAL: sink apparent at node" << endl;
        }
        if ((u>0) && (pSister->u>0) && (pParent->u<0)){ // source apparent
            c_Agcl1 = c_Agcl2 = Cm1(0);
            //             cout << "NON-PHYSICAL: source apparent at node" << endl;
        }

        // define ghost cell values for diffusion
        /*
        if (diffDom){
            c_Dgcl = wSp*c_p1N     + wSs*c_s01;
            c_Dgcr = wSmaj*c_maj01 + wSmin*c_min01;
        }
        else{
            if (u >= 0){
                c_Dgcl = c_Agcl1;
                c_Dgcr = wQmaj*c_maj01 + wQmin*c_min01;
            }
            else{
                c_Dgcl = wQp*c_p1N     + wQs*c_s01;
                c_Dgcr = c_Agcr1;
            }
        }
        */

        c_Dgcl = wSp*c_p1N     + wSs*c_s01;
        c_Dgcr = wSmaj*c_maj01 + wSmin*c_min01;

        // ---------------------------------------------------------------

        if (u >= 0){

            // left b.c.
            rhs(0) += dt*Diff/(h*h)*(1.*c_Dgcl);
            rhs(0) -= dt*u/(2.*h)*(-4.*c_Agcl1 + 1.*c_Agcl2);

            rhs(1) -= dt*u/(2.*h)*(1.*c_Agcl1);

            // right b.c.
            rhs(N) += dt*Diff/(h*h)*(1.*c_Dgcr);
        }

        if (u < 0){

            // left b.c.
            rhs(0) += dt*Diff/(h*h)*(1.*c_Dgcl);

            // right b.c.
            rhs(N) += dt*Diff/(h*h)*(1.*c_Dgcr);
            rhs(N) -= dt*u/(2.*h)*(4.*c_Agcr1 - 1.*c_Agcr2);

            rhs(N-1) -= dt*u/(2.*h)*(-1.*c_Agcr1);
        }
    }

    // Update 'old' concentration value (successive upate)
    //Cm1 = C;

    // solve for concentration:
    if (N > 10){
        C = solverBiCGSTAB.compute((LO).sparseView()).solve(rhs);
    }
    else{
        C = LO.householderQr().solve(rhs);
    }


    // flux at the right border of the last duct
    F_l = Q*C(0) - S*Diff/(2.*h)*(-3.*C(0) + 4.*C(1)   - 1.*C(2));
    F_r = Q*C(N) - S*Diff/(2.*h)*( 3.*C(N) - 4.*C(N-1) + 1.*C(N-2));

    if (reachedEnd){
        cCAW = C(N);
        dcdxCAW = 1./(2.*h)*(3.*C(N) - 4.*C(N-1) + 1.*C(N-2));
    }

    Erel = (LO*C - rhs).norm() / rhs.norm();
    Eabs = (LO*C - rhs).norm();

    if (Erel > 0.1){
        cout << "OOPS: relative error in concentration update is bigger than 10 %" << endl;
    }

    // clean up
    D1_UW *= 0.;
    LO    *= 0.;
    RO    *= 0.;
    rhs   *= 0.;

}


// Update concentration in ducts before going to next time step (global update)
//**************************************************************/
void gas::updateOldConc(){
    Cm1 = C;
}


// Quadratic interpolation with 3 control points and 3 data points (one exact solution)
//**************************************************************/
double gas::quadraticInterpolation(double evalPoint, int opt1, int opt2){

    // declaration
    double c_hat;

    // data points
    if (opt1==0){

        // left side interpolation
        if (opt2==0){
            quadIntC = C.head(3); // with current values
        }
        else{
            quadIntC = Cm1.head(3); // with values from preceeding time-step
        }
    }
    else if (opt1==1){

        // right side interpolation
        if (opt2==0){
            quadIntC = C.tail(3).reverse(); // with current values
        }
        else{
            quadIntC = Cm1.tail(3).reverse(); // with values from preceeding time-step
        }
    }
    else{
      cout << "ERROR: not a valid option for quadratic interpolation of concentration values near bifurcations." << endl;
      return 0.;
    }

    // solve for coefficients
    quadIntA = quadIntM.householderQr().solve(quadIntC);

    // evaluate function
    c_hat = quadIntA(0)*pow(evalPoint,2) + quadIntA(1)*evalPoint + quadIntA(2);

    return c_hat;
}
