function [Po]=Poisson(Dm,lx,ly,m,t,pp,a)
%Poisson Equation Solver using Successive Overrelaxation Method

%Define dimension of the 2-d box 

loop=1;
Nx=size(Dm,1);
Ny=size(Dm,2);

%M=40000; %Maximum iteration value
M=160000;

%Initialize density matrix
Po(1:Nx,1:Ny)=0.0;
des=Dm;

%set boundary Conditions
Po(1,:)=t; 
Po(Nx,:)=t; 
Po(:,1)=t; 
Po(:,Ny)=t;


w=cos(pi/Nx)+cos(pi/Ny); %Converging Term
Ncount=0;
loop=1;

%Choose in which boundaries you would like to set Neumann BC
choice = questdlg('¿In which boundaries will you set Neumann boundary conditions?', ...
	'Neumann boundaries', ...
	'North-South','East-West','East-West');

% Handle response
switch choice
    case 'North-South'
        disp([choice 'The North and South boundaries'])
        
        while loop==1
            Rmin=0; 
            for i=2:Nx-1
                for j=2:Ny-1
                    Residue=w*(0.25.*(Po(i-1,j)+Po(i+1,j)+Po(i,j-1)+Po(i,j+1)+((1/(Nx+1))*(1/(Ny+1)))*des(i,j))-Po(i,j));
                    Rmin=Rmin+abs(Residue);
                    Po(i,j)=Po(i,j)+Residue;
                end
            end
%         
%         u(1,14:42)=u(2,14:42);%NeumannNorth
%         u(Nx,14:42)=u(Nx-1,14:42);%NeumannSouth
            Po(1,:)=Po(2,:);%NeumannNorth
            Po(Nx,:)=Po(Nx-1,:);%NeumannSouth
        
            Rmin=Rmin/(Nx*Ny); % Average Residue per grid point

            ConvC=0.00001; %Convergence term
            %ConvC=0.0001;
            %ConvC=0.001;
            
            if(Rmin>=ConvC) %Check progress
                Ncount=Ncount+1;
                if mod(Ncount,5000) == 0
                    J = 18;
                end
                
                if(Ncount>M)
                    loop=0;
                    disp(['solution doesnt converge in ',num2str(M),' iter'])
                end
            else
                loop=0;
                disp(['solution converges in ',num2str(Ncount) ,' iteration'])
            end
            
        end
        
        
        
        
        
    case 'East-West'
        disp([choice 'East-West boundaries'])

        while loop==1
            Rmin=0; 
            for i=2:Nx-1
                for j=2:Ny-1
                    Residue=w*(0.25.*(Po(i-1,j)+Po(i+1,j)+Po(i,j-1)+Po(i,j+1)+((1/(Nx+1))^2)*des(i,j))-Po(i,j));
                    Rmin=Rmin+abs(Residue);
                    Po(i,j)=Po(i,j)+Residue;
                end
            end

            Po(:,1)=Po(:,2);%NeumannWest
            Po(:,Ny)=Po(:,Ny-1);%NeumannEast

            Rmin=Rmin/(Nx*Ny); % Average Residue per grid point
            ConvC=0.0001;

            if(Rmin>=ConvC)
                Ncount=Ncount+1;
                if(Ncount>M)
                    loop=0;
                    disp(['solution doesnt converge in ',num2str(M),' iter'])
                end
            else
                loop=0;
                disp(['solution converges in ',num2str(Ncount) ,' iteration'])
            end
        end

    case 'East-West'
       disp([choice 'East-West boundaries'])

       while loop==1
           Rmin=0; 
           for i=2:Nx-1
               for j=2:Ny-1
                   Residue=w*(0.25.*(Po(i-1,j)+Po(i+1,j)+Po(i,j-1)+Po(i,j+1)+((1/(Nx+1))^2)*des(i,j))-Po(i,j));
                   Rmin=Rmin+abs(Residue);
                   Po(i,j)=Po(i,j)+Residue;
               end
           end

           Po(:,1)=Po(:,2);%NeumannWest
           Po(:,Ny)=Po(:,Ny-1);%NeumannEast

           Rmin=Rmin/(Nx*Ny); % Average Residue per grid point
           ConvC=0.0001;

           if(Rmin>=ConvC)
               Ncount=Ncount+1;
               if(Ncount>M)
                   loop=0;
                   disp(['solution doesnt converge in ',num2str(M),' iter'])
               end
           else
               loop=0;
               disp(['solution converges in ',num2str(Ncount) ,' iteration'])
           end
       end
        
end %switch

end
