%% The Beginning Part
clc
clear
% Differential Evolution
% Define Control Parameter
Np = 100;
Cr = .5;
F = 0.2;
D = 2;

% for Plotting Function
 x1_p = -10 : 0.1 : 10;
 x2_p = -10 : 0.1 : 10;
 
% set D = 2 for plotting
D_p = 2;
[X1_p,X2_p] = meshgrid(x1_p,x2_p);

% calculate f1_p over the defined meshgrid [X1_p,X2_p] = meshgrid(x1_p,x2_p)
f8_p = X1_p.^2 - 10*cos(2*pi*X1_p) + 10 +...
       X2_p.^2 - 10*cos(2*pi*X2_p) + 10;
% Put initial population into population matirxes
Ip = -10 + (10+10)*rand(D+1,Np);
Pop = Ip; % Permanent Population Matrix
pop_t = zeros(D,Np); % Temporary Population Matrix
    
% Initial Crossover Vector
U = zeros(D,1);
figh = figure;
% MovieVector = zeros(30*D);
    
    for ii = 1 : 30*D
        cla
        
        
        for jj = 1 : Np
            
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
            [f_U,f_X] = RastriginDE(D, U, Pop, jj);
            
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
            
            plot3(Pop(1,jj),Pop(2,jj),Pop(3,jj),'r.','markersize',20);
            
        end
        
        % Put temporary population matrix into permanent population matrix
        Pop(1:D,1:Np) = pop_t;
        best_m = min(Pop(D+1,1:Np));
%         for gg = 1 : Np
%             plot3(Pop(1,gg),Pop(2,gg),Pop(3,gg),'b.','markersize',20);
%             hold on
            
            surfc(X1_p,X2_p,f8_p)

            hold on
            grid on
            view([35 -70])
%             pause(0.1)
            
            MovieVector(ii) = getframe(figh, [10 10 520 400]);
%         end



        
    end
    
      
mywriter = VideoWriter('rastriginde','MPEG-4');
mywriter.FrameRate = 10;
open(mywriter); 
writeVideo(mywriter,MovieVector);
close(mywriter);
    
        


    
    










        
                    
                
            
            
            
            
        
    
    