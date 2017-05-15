figure (1)
%h refine
%force = quadratic, p=1, L21D is direct solver, L21G is gmres
%piecewise2
MeshGrid = [5,10,20,40,80];
E1 = [4.334441e-02,  1.219413e-02, 3.133818e-03, 7.887903e-04, 1.975312e-04];
loglog(MeshGrid,E1,'ro-')
hold on
MeshGrid2 = [5,10,20,40,80];
E2 = [5.275473e-03, 7.023883e-04, 8.981420e-05,1.134002e-05, 1.424493e-06];
loglog(MeshGrid2,E2,'bo-')
hold on
MeshGrid3 = [5,10,20,40,80];
E3 = [3.417899e-04, 2.148093e-05, 1.345928e-06, 8.422031e-08, 5.266802e-09];
loglog(MeshGrid3,E3,'ko-')

E12 = [5.091666e-03, 6.681039e-04, 8.480105e-05, 1.066276e-05, 1.336276e-06];
loglog(MeshGrid3,E12,'b*-');

E13 = [3.437217e-04, 2.153980e-05, 1.347767e-06, 8.427757e-08, 5.268583e-09];
loglog(MeshGrid3,E13,'k*-');
E23 = [3.437217e-04, 2.153980e-05, 1.347767e-06, 8.427757e-08, 5.268583e-09];
loglog(MeshGrid2,E23,'mo-')
legend('p=1,PieceWise2','p=2,PieceWise2','p=3,PieceWise2','VariableP=1,2 PieceWise2','VariableP=1,3 PieceWise2','VariableP=2,3 PieceWise2')
hold on
x = [41,81];
y = [0.4e-3, 0.1e-03];
loglog(x,y,'r');

x = [41,81];
y = [0.000004, 0.5e-6];
loglog(x,y,'b');

x = [41,81];
y = 1.0e-03*[0.00002, 1.25e-06];
loglog(x,y,'k');


%piecewise 3
figure(2)
E1 = [ 7.708352e-02, 2.002275e-02, 5.055198e-03, 1.266941e-03, 3.169325e-04];
loglog(MeshGrid3,E1,'bo-');
hold on
E2=[2.806496e-03, 3.449992e-04, 4.290704e-05, 5.355926e-06, 6.692489e-07];
loglog(MeshGrid3,E2,'ro-');
hold on
E3 = [1.144514e-04, 6.855138e-06, 4.186028e-07, 2.584613e-08, 1.605345e-09];
loglog(MeshGrid3,E3,'ko-');
hold on
%E4 = [1.953592e-15, 2.865145e-15, 7.030268e-15,  1.060437e-14, 1.170625e-13];
%loglog(MeshGrid3,E4,'co-');
hold on
E12 = [3.891643e-02, 1.051095e-02, 2.687032e-03, 6.762493e-04, 1.694240e-04];
loglog(MeshGrid3,E12,'b*-');
hold on
E13 = [3.738724e-02, 1.028231e-02, 2.655324e-03, 6.721116e-04, 1.688989e-04];
loglog(MeshGrid3,E13,'rs-');
hold on
E23 = [1.532467e-03, 1.945891e-04, 2.467534e-05, 3.113747e-06, 3.913428e-07];
loglog(MeshGrid3,E23,'k*-');
hold on
E24 = [1.419892e-03, 1.883488e-04, 2.430869e-05, 3.091564e-06, 3.899795e-07];
loglog(MeshGrid3,E23,'cs-');
hold on
legend('p=1','p=2','p=3','VariableP=1,2','VariableP=1,3','VariableP=2,3','VariableP=2,4')
hold on
x = [41,81];
y = [0.4e-3, 0.1e-03];
loglog(x,y,'b');

x = [41,81];
y = [0.000002, 0.25e-6];
loglog(x,y,'r');

x = [41,81];
y = 1.0e-03*[0.00001, 0.625e-06];
loglog(x,y,'k');