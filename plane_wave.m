% ----- Geometries-----%

x = [-1:.02:1];      % x values in plane 
y = [-1:0.02:1];     % y values in plane
z = 3;      % z values of plane
[X,Y] = meshgrid(x,y);      % Makes plane out of specified x and y points

%-----EM stuff-----%
wv = 0.5;

k = [4, 1, 0];      % k values
K = 2*pi/(wv*(sqrt(k(1)^2+k(2)^2+k(3)^2)))*k;   % Normalised k, so lambda is const

A=0.5;      %Amplitude input
E = zeros((size(X,2)));

for j = 1:size(X,2)
    for l = 1:size(X,2)
         E(j,l) = A*exp(1i*(K(1)*(X(j,l))+K(2)*(Y(j,l))+K(3)*z));  %Creates M with all values of E in there, including phase!
    end
end

% ------ Graphing ----- %
figure(1)
surf(X,Y,angle(E))
shading interp
grid off
colormap('winter')
xlabel('X values')
ylabel('Y Values')
zlabel('Phase')
title('Phase of plane wave in xy plane')
view(3)

figure(2)
surf(X,Y,real(E))
shading interp
grid off
colormap('cool')
xlabel('X values')
ylabel('Y Values')
zlabel('Real Value')
title('Real Part of plane wave in xy plane')
view(3)

figure(3)
surf(X,Y,abs(E))
colormap('spring')
shading interp
grid off
xlabel('X values')
ylabel('Y Values')
zlabel('Real Value')
title('Intensity of plane wave in xy plane')
view(3)