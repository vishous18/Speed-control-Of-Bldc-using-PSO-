clear
clc
warning off;
n =4;             % pop
bird_setp = 4;   % iter
dim = 2;           % variables
c2 = 2;            % PSO parameter C1 
c1 = 0.2;          % PSO parameter C2 
wmax =0.9;         % pso maximum momentum 
wmin =0.4;         % pso minimum momentum
                                               

R1 = rand(dim, n);
% pause
R2 = rand(dim, n);
current_fitness =0*ones(n,1);
rang_Ki=18;
rang_Kp=0.005;
rang_Ki_min=15;
rang_Kp_min=0.0009;

for i=1:dim
    for j=1:n
        if i==1 current_position(i,j)=unifrnd(rang_Ki_min,rang_Ki);
        elseif i==2 current_position(i,j)=unifrnd(rang_Kp_min,rang_Kp);
        end
            fprintf('\t%12.6f ',current_position(i,j));
    end
end

velocity = .3*randn(dim, n) ;
fprintf('\n velocity =%f',velocity);
local_best_position  = current_position ;
                    
for i = 1:n
     Ki=current_position(1,i);
     Kp=current_position(2,i);
     sim('PBLDC');
     xa1=in1.signals.values;
     xa=max(xa1);
     current_fitness(i) = xa;   
     fprintf('\n \t%12.6f \t%12.6f \t%12.6f',Ki,Kp,current_fitness(i));

end

local_best_fitness  = current_fitness; 
[global_best_fitness,g] = min(local_best_fitness) ;
for i=1:n
    global_best_position(:,i) = local_best_position(:,g) ;
end

for iter = 1:bird_setp
    velocity = 0.9 *velocity + c1*(R1.*(local_best_position-current_position)) + c2*(R2.*(global_best_position-current_position));
end
current_position = current_position + velocity;
% current_position
for i=1:dim
    for j=1:n
        if i==1 
            if current_position(i,j)< rang_Ki_min
                current_position(i,j)= rang_Ki_min;
            elseif current_position(i,j)>rang_Ki
                current_position(i,j)=rang_Ki;
            end
        elseif i==2 
            if current_position(i,j)< rang_Kp_min
                current_position(i,j)= rang_Kp_min;
            elseif current_position(i,j)>rang_Kp
                current_position(i,j)=rang_Kp;
            end
        end
    end
end
% current_position

iter = 0 ;        % Iterations’counter
while  ( iter < bird_setp )
iter = iter + 1;
for i = 1:n
     Ki=current_position(1,i);
     Kp=current_position(2,i);
     sim('PBLDC');
      xa1=in1.signals.values;
     xa=max(xa1);
     current_fitness(i) = xa; 
     fprintf(' \n current fit2= %12.6f ',current_fitness(i));
end
for i = 1 : n
        if current_fitness(i) < local_best_fitness(i)
           local_best_fitness(i)  = current_fitness(i);  
           local_best_position(:,i) = current_position(:,i)   ;
        end   
 end
 [current_global_best_fitness,g] = min(local_best_fitness);
    
if current_global_best_fitness < global_best_fitness
   global_best_fitness = current_global_best_fitness;
   for i=1:n
        global_best_position(:,i) = local_best_position(:,g);
   end
end
w = wmax-((wmax-wmin)/bird_setp )*iter;
velocity = w *velocity + c1*(R1.*(local_best_position-current_position)) + c2*(R2.*(global_best_position-current_position));
current_position = current_position + velocity;

% current_position

for i=1:dim
    for j=1:n
        if i==1 
            if current_position(i,j)< rang_Ki_min
                current_position(i,j)= rang_Ki_min;
            elseif current_position(i,j)>rang_Ki
                current_position(i,j)=rang_Ki;
            end
        elseif i==2 
            if current_position(i,j)< rang_Kp_min
                current_position(i,j)= rang_Kp_min;
            elseif current_position(i,j)>rang_Kp
                current_position(i,j)=rang_Kp;
            end
        end
    end
end
fprintf('\n%d  \t%12.6f \t%12.6f \t%12.6f',iter,global_best_fitness,global_best_position(1,1),global_best_position(2,1));
end