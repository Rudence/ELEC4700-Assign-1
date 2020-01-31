% Rudi Hidvary 
% 101037816

close all
clear
clc

nx = 40;
ny = 40;
iterationMax = 50;

V = zeros(ny,nx)
V(:,1) = 1
V(:,end) = 0

% In theory, the current location voltage is calculated as follows for all
% normal cases
% V(i,j) = (1/4)*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1))
% the edge cases can be calculated using a loop method 
% V(:,2) = (1/4)*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1))
for a = 1:iterationMax
    figure(1)
    imagesc(V)
    title('Voltage Plot')
    pause(0.01)
    for i = 1:ny
        for j = 2:nx-1
            % Reinsert the boundary conditions 
            if(i == 1)
                V(i,j) = (1/4)*(2*V(i+1,j)+V(i,j-1)+V(i,j+1))
            elseif(i == ny)
                V(i,j) = (1/4)*(2*V(i-1,j)+V(i,j-1)+V(i,j+1))
            else
                V(i,j) = (1/4)*(V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1))
            end
        end
    end
end 

clc

V = zeros(ny,nx)
V(:,1) = 1
V(:,end) = 1
V(1,:) = 0
V(end,:) = 0 

for a = 1:iterationMax
    figure(2)
    imagesc(V)
    title('Voltage Plot')
    pause(0.01)
    for i = 2:ny-1
        for j = 2:nx-1
            % Reinsert the boundary conditions 
                V(i,j) = (1/4)*(V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1))
        end
    end
end 







            


