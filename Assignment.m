%% The Beginning Part
clc
clear
% Differential Evolution
% Define Control Parameter
Np = 100;
Cr = .5;
F = 0.2;

% User must enter D (Problem's Dimension) as an input
D = input('Please Enter the Problems Dimension: ');
MAX_NFC = 3000 * D;

% Particle Swarm Optimization
% Define Control Parameter
C1 = .01;
C2 = .5;

w = 0.9 : -((0.9-0.4)/(MAX_NFC/Np-1)) : 0.4;
% Initial Population as Ip
%Ip = -10 + (10+10)*rand(D+1,Np);
% Error Matrices
E_DE = zeros (31,MAX_NFC/Np+1);
E_PSO = zeros (31,MAX_NFC/Np+1);

BESTDE = zeros(D+1,31);
BESTPSO = zeros(D+1,31);

for ii = 1 : 31
    
    % Put initial population into population matirxes
    Ip = -10 + (10+10)*rand(D+1,Np);
    Pop = Ip; % Permanent Population Matrix
    pop_t = zeros(D,Np); % Temporary Population Matrix
    
    % Initial Crossover Vector
    U = zeros(D,1);
    NFC = 0; % Reset NFC to Zero at the beginning
    counter = 0;
    

    
    while NFC < MAX_NFC
        counter = counter + 1;
        for jj = 1 : Np
            NFC = NFC + 1;
            parents = randperm(Np,3); % For Xjj, Choose 3 parents randomly 
            
            % Check that none of 3 parents is equal to Xjj
            jjj = [jj,jj,jj];
            while sum(double((parents == jjj))) ~= 0
                parents = randperm(Np,3);
            end
            
            % Mutation
            V = Pop(1:D,parents(1)) + F * (Pop(1:D,parents(2)) - Pop(1:D,parents(3)));
            
            %Crossover
            for kk = 1 : D
                if rand(1) < Cr
                    U(kk,1) = V(kk,1);
                else
                    U(kk,1) = Pop(kk,jj);
                end
            end
            
            %Calculating the f value for X and U vetors
            [f_U,f_X] = WeierstrassDE(D, U, Pop, jj);
            
            Pop(D+1,jj) = f_X; % Insert f_X in the last row of the population matrix
            
            %Choosing next generation of the Xjj vector
            if f_U <= f_X
                X1 = U;
                X11 = 0;
            else
                X1 = Pop(1:D,jj);
                X11 = 1;
            end
            
            % Filling the temporary population matrix 
            pop_t(:,jj) = X1;
            
            % Entering the value of f in the last row 
            if X11 == 0
                Pop(D+1,jj) = f_U;
            else
                Pop(D+1,jj) = f_X;
            end
        end
        
        % Put temporary population matrix into permanent population matrix
        Pop(1:D,1:Np) = pop_t;
        best_m = min(Pop(D+1,1:Np));

    E_DE(ii,counter) = sum(abs(Pop(D+1,:)));    
    end
    [BESTmagnitude,BESTindex] = min(Pop(D+1,1:Np));
    BESTDE(:,ii) = Pop(:,BESTindex);
    E_DE(ii,MAX_NFC/Np+1) = sum(E_DE(ii,:));
    % Set Initial Matrixes and Parameters
    Pop = Ip; % Initial Population Matrix
    Pbest = Pop; % Initial Pbest Matrix
    v = zeros(D,Np); % Initial V Matirx
    NFC = 0; % Reset NFC to Zero at the beginning
    w_index = 0;
    counter = 0;
    for jj = 1 : Np
        f_X = WeierstrassPSO(D, Pop, jj);
        Pop(D+1,jj) = f_X;
        Pbest(D+1,jj) = f_X;
    end
    
    % Finding Gbest for the Initial Population 
    [MINmagnitude,MINindex] = min(Pbest(D+1,1:Np));
    Gbest = Pbest(:,MINindex);
    
        
    while NFC < MAX_NFC
        w_index = w_index + 1;
        counter = counter + 1;
        for kk = 1 : Np
            NFC = NFC + 1;
            va = w(w_index) * v(1:D,kk) + C1 * (Pbest(1:D,kk)-Pop(1:D,kk)) + C2 * (Gbest(1:D)-Pop(1:D,kk));
            
            for ll = 1 : D
                if va(ll) < -0.2
                    va(ll) = -0.2;
                elseif va(ll) > 0.2
                    va(ll) = 0.2;
                end
            end
                      
            Pop(1:D,kk) = Pop(1:D,kk) + va;
            v(1:D,kk) = va;
            
            
                        
            f_X = WeierstrassPSO(D, Pop, kk);
            Pop(D+1,kk) = f_X; % Insert f_X in the last row of the population matrix
            
            if f_X <= Pbest(D+1,kk)
                Pbest(1:D,kk) = Pop(1:D,kk);
                Pbest(D+1,kk) = f_X;
            end
            if f_X <= Gbest(D+1)
                Gbest(1:D) = Pop(1:D,kk);
                Gbest(D+1) = f_X;
                Gbest(D+1);
                
            end
        end

    E_PSO(ii,counter) = sum(abs(Pop(D+1,:)));
    end
    [BESTmagnitude,BESTindex] = min(Pop(D+1,1:Np));
    BESTPSO(:,ii) = Pop(:,BESTindex);
    E_PSO(ii,MAX_NFC/Np+1) = sum(E_PSO(ii,:));
    
    
end

[MINEPSOmagnitude,MINEPSOindex] = min(E_PSO(:,MAX_NFC/Np+1));
MINEPSOmagnitude

[MINEDEmagnitude,MINEDEindex] = min(E_DE(:,MAX_NFC/Np+1));
MINEDEmagnitude

[MAXEPSOmagnitude,MAXEPSOindex] = max(E_PSO(:,MAX_NFC/Np+1));
MAXEPSOmagnitude
[MAXEDEmagnitude,MAXEDEindex] = max(E_DE(:,MAX_NFC/Np+1));
MAXEDEmagnitude

meanPSO = (sum(E_PSO(:,MAX_NFC/Np+1)))/31
meanDE = (sum(E_DE(:,MAX_NFC/Np+1)))/31

STDPSO = std(E_PSO(:,MAX_NFC/Np+1))
STDDE = std(E_DE(:,MAX_NFC/Np+1))

X = Np : Np : MAX_NFC;
YDE = E_DE(MINEDEindex, 1:MAX_NFC/Np);
YPSO = E_PSO(MINEPSOindex, 1:MAX_NFC/Np);
p = plot(X,YDE,'b--',X,YPSO);
p(1).LineWidth = 2;
p(2).LineWidth = 2;
grid on
legend('DE','PSO')
title(['Best Fitness Error over 31 Runs, D = ',num2str(D)])

[THEBESTDEmagnitude,THEBESTDEindex] = min(BESTDE(D+1,1:31));
   
BESTDE(:,THEBESTDEindex)

[THEBESTPSOmagnitude,THEBESTPSOindex] = min(BESTPSO(D+1,1:31));
   
BESTPSO(:,THEBESTPSOindex)



        
                    
                
            
            
            
            
        
    
    