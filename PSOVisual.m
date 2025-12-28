%% The Beginning Part
clc
clear
% Differential Evolution
% Define Control Parameter
D = 2;
Np = 100;
MAX_NFC = 3000 * D;

% User must enter D (Problem's Dimension) as an input


% for Plotting Function
 x1_p = -10 : 0.1 : 10;
 x2_p = -10 : 0.1 : 10;
 
% set D = 2 for plotting
D_p = 2;
[X1_p,X2_p] = meshgrid(x1_p,x2_p);

% calculate f1_p over the defined meshgrid [X1_p,X2_p] = meshgrid(x1_p,x2_p)
f8_p = X1_p.^2 - 10*cos(2*pi*X1_p) + 10 +...
       X2_p.^2 - 10*cos(2*pi*X2_p) + 10;

% Particle Swarm Optimization
% Define Control Parameter
C1 = .01;
C2 = .5;

w = 0.9 : -((0.9-0.4)/(MAX_NFC/(Np)-1)) : 0.4;
% Initial Population as Ip
%Ip = -10 + (10+10)*rand(D+1,Np);
% Error Matrices





    
    % Put initial population into population matirxes
    Ip = -10 + (10+10)*rand(D+1,Np);
    % Set Initial Matrixes and Parameters
    Pop = Ip; % Initial Population Matrix
    Pbest = Pop; % Initial Pbest Matrix
    v = zeros(D,Np); % Initial V Matirx
    w_index = 0;

     
    for jj = 1 : Np
        f_X = RastriginPSO(D, Pop, jj);
        Pop(D+1,jj) = f_X;
        Pbest(D+1,jj) = f_X;
    end
    
    % Finding Gbest for the Initial Population 
    [MINmagnitude,MINindex] = min(Pbest(D+1,1:Np));
    Gbest = Pbest(:,MINindex);
figh = figure;    
for ii = 1 : 30*D
    cla 
    w_index = w_index + 1;

        for kk = 1 : Np
            
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
            
            
                        
            f_X = RastriginPSO(D, Pop, kk);
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
            plot3(Pop(1,kk),Pop(2,kk),Pop(3,kk),'r.','markersize',20);

        end
        surfc(X1_p,X2_p,f8_p)

        hold on
        grid on
        view([-120 70])
%             pause(0.1)
            
        MovieVector(ii) = getframe(figh, [10 10 520 400]);

   
end
    
mywriter = VideoWriter('rastrigin','MPEG-4');
mywriter.FrameRate = 10;
open(mywriter); 
writeVideo(mywriter,MovieVector);
close(mywriter);    










        
                    
                
            
            
            
            
        
    
    