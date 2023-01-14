
N = 100;
T = linspace(300, 600, 100);
P = linspace(1, 60, 100);
equation = 'gc';
molecule = 'benzene';
Hr=zeros(N,N);
Sr = zeros(N,N);
Gr = zeros(N,N);

for i=1:N
    for j=1:N
        state.T = T(j);
        state.P = P(i);
        [Hr(i,j), Sr(i,j)] = residual(state,molecule,equation);
        Gr(i,j) = Hr(i,j) - state.T*Sr(i,j);
    end
end
    