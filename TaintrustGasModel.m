function NewGasModel()
 %Code for a network of gas pipelines with compressors. 
% The code includes multiple functions for calculating gas flow through pipes and compressors, adjusting pressures, 
% and solving the system iteratively to reach aerging solution

close  %close closes any open figure windows.
clc    %clears the command window for a clean start.

%% DATA
Model_Compressor=1; %flag to enable or disable compressor, 1 is enable 0 is disable

[Branch_Connectivity,Branch_Length,Branch_Diameter,Branch_Efficiency,...  %get gas predefined data related to the gas network into function GetGasNetworkData  
    Node_Demand,Compressor_Location,Compressor_Efficiency,...
    Compressor_Costs,Compressor_Ratio,Standard_Temp,Standard_Pressure,...
    Average_Temperature,Gas_Gravity,Gas_Comp,Compressor_Temp,...
    Gas_Heat,Compressor_iComp]=GetGasNetworkData();
Slack=1;                                        %Slack is a reference node with a fixed pressure (slack node), used to stabilize the network.
SlackP=1000;%psia                    %SlackP is the pressure at the slack node, set to 1000 psia.

Pressure=[%psia  Initial pressure estimates for each node in the system, in psia (pounds per square inch absolute)
    1000   %There are 14 nodes in total  node 1 pressure here is a placeholder
    500
    600
    400
    500
    450
    300
    670
    900
    800
    700
    500
    600
    550
    ];
% Branch_Connectivity=Branch_Connectivity(1:5,:);Branch_Connectivity(5,:)=[3 5]
% Compressor_Location=[4;5];
% Compressor_Ratio=Compressor_Ratio(1:2,:);
% Node_Demand=zeros(5,1);Node_Demand(4:5)=[215;148]; 



%% Calculations
Ng=max(max(Branch_Connectivity));  %Ng is total node
Nbg=length(Branch_Connectivity(:,1));   %No. of branches (Pipes)
Ncg=length(Compressor_Location(:,1));   %Total compressor  where,  (:) specify all the rows of the first column (1)



%%The purpose of these lines is to assign specific values into the matrix Node_aux_Full, based on data stored in other matrices
% or arrays, specifically Branch_Connectivity, Compressor_Location, and Compressor_Ratio.
Node_aux_Full=zeros(Ng,5);   %Node_aux_Full will be a matrix with Ng rows and 5 columns, filled with zeros. 

 %extracts specific elements from the first column of Branch_Connectivity based on the Compressor_Location. 
 %this line is populating the first column of Node_aux_Full with indices from 1:Ncg in specific rows, 
 % determined by the connectivity information in Branch_Connectivity.
Node_aux_Full(Branch_Connectivity(Compressor_Location,1),1)=1:Ncg;

%Location of additional demand to be added to those nodes
%This line assigns the values from the second column of Branch_Connectivity (at rows specified by Compressor_Location) to the second 
% column of Node_aux_Full at the same rows specified by Branch_Connectivity(Compressor_Location,1).
Node_aux_Full(Branch_Connectivity(Compressor_Location,1),2)=... %Fill column 2
    Branch_Connectivity(Compressor_Location,2);
%Nodes with pressure affected by compressors
Node_aux_Full(Branch_Connectivity(Compressor_Location,2),3:4)=...
    Compressor_Ratio;
%Node to be considered;
Node_aux_Full(Branch_Connectivity(Compressor_Location(...
    Compressor_Ratio(:,2)~=0),2),5)=Branch_Connectivity(...
    Compressor_Location(Compressor_Ratio(:,2)~=0),1);
%Compressor_Ratio(:,2)~=0 creates a logical index, selecting rows where the second column of Compressor_Ratio is not zero.
Pressure(Slack)=SlackP;  %node 1 is slack node
ConvFlg=Inf;
xcou=0;

%main program


while ConvFlg > 0.0001 %keep repeating until error is less than 0.0001
    xcou=xcou+1;
    %Get flows at each node and thorugh each pipe
    %[New_Branch_Flow,New_Node_Flow,New_CompressorIn,Pipe_Location] are outputs of GetGasFlows function

    Branch_Diameter([1, 2]) = [22, 21]; % Increase diameters for Branch 1-2 and 1-3 in inches
    Branch_Length([1, 2]) = [90, 95]; % Increase pipeline lengths for Branch 1-2 and 1-3 in miles

    [New_Branch_Flow,New_Node_Flow,New_CompressorIn,Pipe_Location]=GetGasFlows(Pressure,...
        Ng,Nbg,Ncg,Compressor_Location,Branch_Length,Branch_Diameter,...
        Branch_Connectivity,Branch_Efficiency,Standard_Temp,...
        Standard_Pressure,Gas_Gravity,Average_Temperature,...
        Gas_Comp,Compressor_Temp,Compressor_Efficiency,Gas_Heat,...
        Compressor_iComp,Compressor_Costs,Model_Compressor,...
        Compressor_Ratio,Node_aux_Full);

    
      % Power flow in 1-2 and 1-3
% Combine gas flow from pipe 1-2 and pipe 1-3 (es 1 and 2)
flow_12 = abs(New_Branch_Flow(1)); % Gas flow in pipe 1-2
flow_13 = abs(New_Branch_Flow(2)); % Gas flow in pipe 1-3
flow_23 = abs(New_Branch_Flow(3));
flow_36 = abs(New_Branch_Flow(6));


% Combine the gas flows
combined_flow_MMcfd = flow_12 + flow_13; % Combined flow in million cubic feet/day
combined_consump = flow_13 + flow_23 - flow_36;    %Gas node 3 consumption
% Convert to cubic meters per second
combined_flow_m3ps = combined_flow_MMcfd * 1e6 * 0.0283168 / (86400) ;
combined_flow_node3 = combined_consump * 1e6 * 0.0283168 / (86400) ;

% Use GCV to calculate power flow
GCV = 39.5/0.87; % MJ/cubic meter
power_flow_MJ = combined_flow_m3ps * GCV/1000;
power_flow_node3 = combined_flow_node3 * GCV/0.9/1000; %in actual should not include /1000

% Display the results
fprintf('Combined gas flow in pipes 1-2 and 1-3: %.4f million cubic feet/day\n', combined_flow_MMcfd);
fprintf('Combined gas flow in cubic meters per second: %.6f m^3/s\n', combined_flow_m3ps);
fprintf('Power flow P2G (P = q * GCV): %.2f MW\n\n', power_flow_MJ);
fprintf('Combined gas flow in cubic meters per second of node 3: %.4f m^3/s\n', combined_flow_node3);
fprintf('Power flow of Generator (node 3) (P = q * GCV): %.2f MW\n\n', power_flow_node3);



    %Get difference between expected and simulated nodal flows
    aux=Branch_Connectivity(Compressor_Location,:);
    aux_Demand=Node_Demand; %aux_Demand: This is ified version of Node_Demand, the gas demand at each node.
    aux_Demand(aux(:,1)) = aux_Demand(aux(:,1)) + aux_Demand(aux(:,2));
    aux_Demand(aux(:,2)) = 0;  %The demand at the downstream node is set to zero because the compressor will meet that demand, effectively shifting it upstream.

    aux=Branch_Connectivity(Compressor_Location,:);
    New_Node_Flow(aux(:,1)) = New_Node_Flow(aux(:,1)) + New_Node_Flow(aux(:,2));
    New_Node_Flow(aux(:,2)) = 0;
    New_DF = New_Node_Flow+aux_Demand; %includes both the original node flows and the compressor-induced flows.
    New_DF(aux(:,1))=New_DF(aux(:,1)) + New_CompressorIn;
    
    
    Real_Nodes = 1:Ng;
    Real_Nodes = Real_Nodes(Node_aux_Full(:,4)==0);   %excludes nodes that are not affected by compressors
    
    %Removes the slack node from the list of nodes whose pressures will be updated. 
    % The slack node has a fixed pressure, so it doesn't need to be included in the pressure update calculations.
    srn=length(Real_Nodes);
    aux=find(Real_Nodes==Slack);
    if aux==1 %If the slack node is the first node, then Real_Nodes_NoS will be all the nodes except the first one.
        Real_Nodes_NoS=Real_Nodes(2:srn);
    elseif Slack==srn  %If the slack node is the last node, then Real_Nodes_NoS will exclude the last node.
        Real_Nodes_NoS=Real_Nodes(1:srn-1);
    else  %Otherwise, Real_Nodes_NoS excludes the slack node from the list of real nodes.
        Real_Nodes_NoS=Real_Nodes([1:aux-1 aux+1:srn]);
    end

    %Build Jacobian
    Jg=GetGasJacobian(Ng,Slack,Pipe_Location,Branch_Connectivity,Pressure,...
        Branch_Efficiency,Branch_Length,Branch_Diameter,Standard_Temp,...
        Standard_Pressure,Gas_Gravity,Average_Temperature,Gas_Comp,Ncg,...
        Compressor_Location,New_Branch_Flow,New_CompressorIn,...
        Compressor_iComp,Gas_Heat,Compressor_Temp,Compressor_Efficiency,...
        Compressor_Costs,Compressor_Ratio,Model_Compressor,Real_Nodes,...
        Node_aux_Full);

    
    %Calculate pressure updates
    %Iterative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Pressure_Correction=Jg\New_DF(Real_Nodes_NoS);  % "\" is matrix division
    %New_DF represents the flow mismatch (how much the flow at each node differs from what is expected
    %Update convergence criteria
    ConvFlg=max(abs(Pressure_Correction)); %When ConvFlg becomes smaller than 0.0001, the system is considered to have converged, meaning the pressure values are stable.
    Pressure(Real_Nodes_NoS)= Pressure(Real_Nodes_NoS) + Pressure_Correction;
    %This line adjusts the pressures based on the corrections calculated in the previous step, moving the system closer to the solution where flows and pressures are balanced.

    if xcou>200
        error('Focing termination'); %If the loop hits 100 iterations without converging, it throws an error and terminates the simulation.
    end
end

%     Display pressures
fprintf('Recommended pressures at each node:\n');
for x1=1:Ng
    fprintf('P:%2.0f:  %f\n',x1,Pressure(x1))
end
fprintf('Maximum error: %.4f\n',ConvFlg);
fprintf('Iterations: %.d\n',xcou);
% tstaux=input('');


% Display compressor flows, connections, and gas consumption
fprintf('\nCompressor gas flows, connections, and consumption:\n');
for xc = 1:Ncg
    % Get nodes connected by each compressor
    node_from = Branch_Connectivity(Compressor_Location(xc), 1);
    node_to = Branch_Connectivity(Compressor_Location(xc), 2);
    
    % Retrieve the flow rate through the compressor from New_Branch_Flow
    compressor_flow_rate = abs(New_Branch_Flow(Compressor_Location(xc)));  % Flow in million cubic feet per day

    % Get the gas consumption from New_CompressorIn
    gas_consumption = New_CompressorIn(xc);  % already in million cubic feet per day

    % Output compressor details
    fprintf('Compressor %d (Node %d -> Node %d):\n', xc, node_from, node_to);
    fprintf('  Flow Rate through Compressor = %.4f million cubic feet per day\n', compressor_flow_rate);
    fprintf('  Gas Consumption = %.4f million cubic feet per day\n', gas_consumption);
end


%% Get gas flows
function [New_Branch_Flow,New_Node_Flow,New_CompressorIn,Pipe_Location]=GetGasFlows(...
    Pressure,Ng,Nbg,Ncg,Compressor_Location,Branch_Length,...
    Branch_Diameter,Branch_Connectivity,Branch_Efficiency,Standard_Temp,...
    Standard_Pressure,Gas_Gravity,Average_Temperature,Gas_Comp,...
    Compressor_Temp,Compressor_Efficiency,Gas_Heat,Compressor_iComp,...
    Compressor_Costs,Model_Compressor,Compressor_Ratio,Node_aux_Full)

ax1=ones(Nbg,1);ax1(Compressor_Location)=0;
ax2=1:Nbg;
Pipe_Location=ax2(ax1==1);

%Calculate flow throughout pipes
New_Branch_Flow=zeros(Nbg,1);
New_Node_Flow=zeros(Ng,1);
for x1=Pipe_Location
    if x1==1
        clc
        fprintf('There are 12 pipelines, 14 nodes and 4 compressors\n')
        fprintf('Gas flow in each branch:\n')
    end
    L=Branch_Length(x1);  %Miles
    D=Branch_Diameter(x1);  %inches
    
    %Friction factor
    F=0.128/D^(1/3);
    
    
    %Get pressures
    Bc = Branch_Connectivity(x1,:);
    p=[0 0];
    for x2=1:2
        switch Node_aux_Full(Bc(x2),4)
            case 0
                p(x2) = Pressure((Bc(x2)));
            case 1
                p(x2) = Pressure(Node_aux_Full(Bc(x2),5)) * Node_aux_Full(Bc(x2),3);
            case 2
                p(x2) = Node_aux_Full(Bc(x2),3);
        end
    end
    
%     p1=Pressure(Branch_Connectivity(x1,1));%psia
%     p2=Pressure(Branch_Connectivity(x1,2));%psia
    
    Pipe_Eff=Branch_Efficiency(x1);
%     Pipe_Eff=1;
    
    %Get flow throughout the pipe
   ax= SignP(p(1),p(2)) * Pipe_Eff * 6.4774 * Standard_Temp * D^2.5 / ...
        Standard_Pressure / sqrt(F*Gas_Gravity*L*Average_Temperature*...
        Gas_Comp) * sqrt( SignP(p(1),p(2))*(p(1)^2-p(2)^2) );
%     ax= SignP(p(1),p(2)) * Pipe_Eff * 6.4774 * Average_Temperature * D^2.5 / ...
%         Standard_Pressure  / sqrt(F*Gas_Gravity*L*Standard_Temp*...
%         Gas_Comp) * sqrt( SignP(p(1),p(2))*(p(1)^2-p(2)^2) );
    New_Branch_Flow(x1)=ax*24/1000000;%ax in cubic feet per hour and new_branch_flow is in million cubic feet per day
    
    fprintf('(Node %2.0f -> Node%2.0f)%10.4f - %10.4f: %10.4f\n',Bc, p(1),p(2),New_Branch_Flow(x1))    %New_Branch_Flow(x1)
%     if x1 ==1
%         Pipe_Eff * 6.4774 * Standard_Temp * D^2.5 / ...
%         Standard_Pressure / sqrt(F*Gas_Gravity*L*Average_Temperature*...
%         Gas_Comp) * sqrt( SignP(p(1),p(2))*(p(1)^2-p(2)^2) )*24/1000000
%         
%         0.9 * 6.4774 * 492.3000 * 19.6^2.5 / ...
%         14.5038 / sqrt(0.0475*0.6*80.50*520*...
%         0.87) * sqrt( (1000^2-676.65^2) )*24/1000000
% 
%         p(2)
%     end

    %Estimate nodal flows
    New_Node_Flow(Branch_Connectivity(x1,1)) = ...
        New_Node_Flow(Branch_Connectivity(x1,1)) + New_Branch_Flow(x1);
    New_Node_Flow(Branch_Connectivity(x1,2)) = ...
        New_Node_Flow(Branch_Connectivity(x1,2)) - New_Branch_Flow(x1);
end
fprintf('\n');
%Calculate flows throughout compressors
New_CompressorIn=zeros(Ncg,1);
if Model_Compressor==1
    xc=0;
    for x1 = Compressor_Location'
        xc=xc+1;
        p1=Pressure(Branch_Connectivity(x1,1));%psia
        if Compressor_Ratio(xc,2)==1
            p2=p1*Compressor_Ratio(xc,1);
        else
            p2=Compressor_Ratio(xc,1);
        end

        B = 0.08530992*Compressor_Temp/Compressor_Efficiency(xc)*(Gas_Heat/(Gas_Heat-1));
        
        f=-New_Node_Flow(Branch_Connectivity(x1,1))/24*1000000;%
%********* calculates the horsepower (HP) required by the compressor to increase the gas pressure from p1 to p2
        HP = B * f * ( (p2/p1)^(Compressor_iComp*((Gas_Heat-1)/Gas_Heat)) - 1);
        
        Ca=Compressor_Costs(xc,1);
        Cb=Compressor_Costs(xc,2);
        Cc=Compressor_Costs(xc,3);
%calculates the amount of gas "injected" by the compressor into the system, based on the operating cost and the horsepower required to compress the gas
        New_CompressorIn(xc)=(Ca+Cb*HP+Cc*HP^2)*24/1000000;%ft^3/h --> million cubic feet per day Gas consumption rate
% gas flow in branch after the compressor is adjusted by subtracting the compressor's contribution (New_CompressorIn) from the original flow rate (f).        
        New_Branch_Flow(x1) = f/1000000*24 - New_CompressorIn(xc);  %Compressor flow rate
  
        %Estimate nodal flows
        New_Node_Flow(Branch_Connectivity(x1,1)) = ...
            New_Node_Flow(Branch_Connectivity(x1,1)) + New_Branch_Flow(x1); %flow at the upstream node (node 1) is increased by the adjusted branch flow
        New_Node_Flow(Branch_Connectivity(x1,2)) = ...
            New_Node_Flow(Branch_Connectivity(x1,2)) - New_Branch_Flow(x1);
    end
end



%% Get gas Jacobian   allowing the simulation to adjust pressures at each node iteratively until the system converges.
function Jg=GetGasJacobian(Ng,Slack,Pipe_Location,Branch_Connectivity,...  
    Pressure,Branch_Efficiency,Branch_Length,Branch_Diameter,...
    Standard_Temp,Standard_Pressure,Gas_Gravity,Average_Temperature,...
    Gas_Comp,Ncg,Compressor_Location,New_Branch_Flow,New_CompressorIn,...
    Compressor_iComp,Gas_Heat,Compressor_Temp,Compressor_Efficiency,...
    Compressor_Costs,Compressor_Ratio,Model_Compressor,Real_Nodes,...
    Node_aux_Full)

Ng=length(Real_Nodes);
x1=1;
while Slack~=Real_Nodes(x1)  %The code will stop when find slack node
    x1=x1+1;
end
if x1==1  %Here, removing the slack node from the list of nodes
    Xng=Real_Nodes(2:Ng);
elseif x1==Ng    % The actual code say x2 == here
    Xng=Real_Nodes(1:Ng-1);
else
    Xng=Real_Nodes([1:Slack-1 Slack+1:Ng]);
end

Jg=zeros(Ng-1,Ng-1);  % initializes the Jacobian matrix (Jg)
Dp=0.000001;   %Dp is a small value used for numerical calculations.


%Auxiliar for assigning values to the Jacobian
xJaux=zeros(max(Xng),1);  
xJaux(Xng)=1:Ng-1;

%Get Jacobian elements corresponding to flow through the pipes
for xp = Pipe_Location%For each pipe
    
    %Get pipe efficiency
    Pipe_Eff=Branch_Efficiency(xp);
    
    %Get pipe length
    L=Branch_Length(xp);%Miles
    
    %Get pipe internal diameter
    D=Branch_Diameter(xp);%inches
    
    %Get friction factor
    F=0.128/D^(1/3);
    
    %Identify relevant nodes and get pressures
    Bc = Branch_Connectivity(xp,:);
    
    p=[0 0];
    Maux = [1 1];   %Use to adjust the pressure
    Eaux=[1 1];
    for x1=1:2
        switch Node_aux_Full(Bc(x1),4)
            %The pressure is based on this node
            case 0
                %Get pressure
                p(x1) = Pressure((Bc(x1)));
                
                
                %Pressure is a function of a different node
            case 1
                %Add auxiliar
                Maux(x1) = Node_aux_Full(Bc(x1),3);
                
                %Change node location
                Bc(x1)=Node_aux_Full(Bc(x1),5);
                
                %Get pressure
                p(x1) = Pressure(Bc(x1)) * Maux(x1);
                
                
                %Is pressure fixed at this node?
            case 2
                %Get pressure
                p(x1) = Node_aux_Full(Bc(x1),3);
                
                %Change node location
                Bc(x1)=Node_aux_Full(Bc(x1),5);
                
                %Do not differentiate this side
                Eaux(x1)=0;
        end
    end
    
    
%     xp
%     Bc
%     error('Just stop');

    %Differentiate input and output sides
    aux = [1 -1];
    ax1 = Pipe_Eff * 6.4774 * Standard_Temp * D^2.5 / Standard_Pressure ...
        / sqrt( F * Gas_Gravity * L * Average_Temperature * Gas_Comp) ...
        *24/1000000;
    ax2 = sqrt( SignP(p(1),p(2)) * (p(1)^2 - p(2)^2) );
    for x1=1:2
%         fprintf('First ');
%         Node_aux_Full(Bc(1),4)
        %Differentiate if the pressure is not fixed and it is not the slack
        if Node_aux_Full(Bc(x1),4)~=2 && xJaux(Bc(x1))~=0 && Eaux(x1)==1
%             fprintf(' Second');
            Jval = aux(x1) * ax1/ax2 * Maux(x1) * p(x1) ;
            
            %Store results
            ax3=[-1 1];
            for x2=1:2
                if xJaux(Bc(x2))~=0
% fprintf('Df%d_%dp%d: %f --> J(%d,%d)\n',Bc,Bc(x1),Jval* ax3(x2),xJaux(Bc(x2)),xJaux(Bc(x1)))
% error('Just stop');
                    Jg(xJaux(Bc(x2)),xJaux(Bc(x1))) = ...
                        Jg(xJaux(Bc(x2)),xJaux(Bc(x1))) + Jval * ax3(x2);
                end
            end
        end
%         fprintf('\n');
    end
%     Jg
%     error('Just stop');
end



%Get Jacobian elements corresponding to compressors
if Model_Compressor==1
    for xc=1:Ncg  %Ncg == 4
%For each compressor xc, Bc retrieves the nodes connected by the compressor. The switch statement manages different compressor behaviors depending on Compressor_Ratio(xc, 2).
        Bc = [
            Branch_Connectivity(Compressor_Location(xc),1)
            Branch_Connectivity(Compressor_Location(xc),2)
            ];
        switch Compressor_Ratio(xc,2)
            %Is the pressure a function of p1
            case 1
                aux = Compressor_Ratio(xc,1);
                Jval = ( 0.08530992 * Compressor_Temp * (Gas_Heat/(Gas_Heat-1))...
                    / Compressor_Efficiency(xc) * ( aux^ ( Compressor_iComp * ...
                    (Gas_Heat-1) / Gas_Heat )-1 ) ) * Compressor_Costs(xc,2) * ...
                    -sum(Jg(1:xJaux(Bc(1))-1,xJaux(Bc(1))));
                
            case 2
                %Use numeric integration for the time being
                %Look for pipes connected to this node
                f=[0 0];Dp=[0 0.00001];
                Ps=Pressure;
                for xf=1:2
                    Ps(Bc(1))=Ps(Bc(1))+Dp(xf);
                    for xp = Pipe_Location%For each pipe
                        %Identify relevant nodes
                        Bcp = Branch_Connectivity(xp,:);
                        if Bcp(1)==Bc(1) || Bcp(2)==Bc(1)
                            p=[0 0];
                            for x1=1:2
                                switch Node_aux_Full(Bcp(x1),4)
                                    case 0
                                        p(x1) = Ps((Bcp(x1)));
                                    case 1
                                        p(x1) = Ps(Node_aux_Full(Bcp(x1),5)) * Node_aux_Full(Bcp(x1),3);
                                    case 2
                                        p(x1) = Node_aux_Full(Bcp(x1),3);
                                end
                            end
                            Pipe_Eff=Branch_Efficiency(xp);
                            L=Branch_Length(xp);%Miles
                            D=Branch_Diameter(xp);%inches
                            F=0.128/D^(1/3);
                            
                            
                            aux = SignP(p(1),p(2)) * Pipe_Eff * 6.4774 * ...
                                Standard_Temp * D^2.5 / Standard_Pressure / ...
                                sqrt( F * Gas_Gravity * L * ...
                                Average_Temperature * Gas_Comp) *24/1000000*...
                                sqrt( SignP(p(1),p(2)) * (p(1)^2 - p(2)^2) );
                            
                            if Bcp(1)==Bc(2)
                                f(xf)=f(xf)+aux;
                            else
                                f(xf)=f(xf)-aux;
                            end
                        end
                    end
                end                                

                p=[Pressure(Bc(1)) Compressor_Ratio(xc,1)];
                f(1)= 0.08530992 * Compressor_Temp * (Gas_Heat/(Gas_Heat-1))... %flow through the compressor at the current upstream pressure p(1).
                    / Compressor_Efficiency(xc) * ( (p(2)/p(1))^ ( Compressor_iComp * ...
                    (Gas_Heat-1) / Gas_Heat )-1 ) * Compressor_Costs(xc,2) * f(1);
                
                p=[Pressure(Bc(1))+Dp(2) Compressor_Ratio(xc,1)];
                f(2)= 0.08530992 * Compressor_Temp * (Gas_Heat/(Gas_Heat-1))...
                    / Compressor_Efficiency(xc) * ( (p(2)/p(1))^ ( Compressor_iComp * ...
                    (Gas_Heat-1) / Gas_Heat )-1 ) * Compressor_Costs(xc,2) * f(2);
                Jval = -(f(2)-f(1))/Dp(2);   
        end
        

% fprintf('Df%d_%dp%d: %f --> J(%d,%d)\n',Bc,Bc(1),Jval* ax3(x2),xJaux(Bc(1)),xJaux(Bc(1)))        
%         fprintf('%d-%d/%d: %f \n',Bc,Bc(1),Jval)
%         error('Just stop');
        
        %Add results to the Jacobian
        Jg(xJaux(Bc(1)),xJaux(Bc(1))) = ...
            Jg(xJaux(Bc(1)),xJaux(Bc(1))) - Jval;
    end
end
%  error('Just stop');


%% Colebrook equation  calculates the friction factor based on the Reynolds number and pipe roughness
function f=Colebrook(err,Re)
f=0.01;
dt=0.00000001;

faux=[Inf f]; %Store previous and current values of f for iteration
x1=0;
while abs(faux(1)-faux(2))>dt  %function iterates until the change in f (the friction factor) is less than a small tolerance dt
    x1=x1+1;
    %Assessing the function
    ax1=-2*log(err/3.7+2.51/Re/sqrt(f))-1/sqrt(f);
        
    %Numeric differentiation
    ax2=-2*log(err/3.7+2.51/Re/sqrt(f+dt))-1/sqrt(f+dt);
    ax2=(ax2-ax1)/dt;
    
    %Avoid using negative values
    f=max(f-ax1/ax2,dt);
    
    faux=[faux(2) f];
    
    if x1==100
        fprintf('A solution for the Colebrook equation could not \n');
        fprintf('be found after %d itertions\n',x1);
        fprintf('err:%.4f; Re:%.4f\n',err,Re)
        error('Stopping program');
    end
end

%% Sign function
function no=SignP(a,b)

if a>=b
    no=1;
else
    no=-1;
end

%% Network data
function [Branch_Connectivity,Branch_Length,Branch_Diameter,Branch_Efficiency,...
    Node_Demand,Compressor_Location,Compressor_Efficiency,...
    Compressor_Costs,Compressor_Ratio,Standard_Temp,Standard_Pressure,...
    Average_Temperature,Gas_Gravity,Gas_Comp,Compressor_Temp,...
    Gas_Heat,Compressor_iComp]=GetGasNetworkData()

Branch_Connectivity=[
    1 2   %1
    1 3   %2
    2 3   %3
    2 4   %4
    4 5   %5
    3 6   %6
    6 7   %7
    5 8   %8
    8 9   %9
    7 10  %10
    10 11 %11
    9 12  %12
    11 13 %13
    12 13 %14
    12 14 %15
    13 14 %16
    ];

Branch_Length=[%Miles
    80.5 %1
    88.3 %2
    55.9 %3
    61.1 %4
    0    %5
    67.9 %6
    0    %7
    93.5 %8
    0    %9
    79.7 %10
    0    %11
    73.5 %12
    87.9 %13
    86.6 %14
    79.7 %15
    83.5 %16
    ];
Branch_Diameter=[%inch
    19.6 %1 
    19.6 %2
    19.6 %3
    19.6 %4
    0    %5 0
    19.6 %6
    0    %7 0
    19.6 %8
    0   %9 0
    16.7 %10
    0    %11 0
    16.7 %12
    16.7 %13
    16.7 %14
    16.7 %15
    16.7 %16
    ];

Branch_Efficiency=[
    0.9  %1
    0.9  %2
    0.9  %3
    0.9  %4
    0.9  %5
    0.9  %6
    0.9  %7
    0.85 %8
    0.85 %9
    0.9  %10
    0.9  %11
    0.85 %12
    0.85 %13
    0.9  %14
    0.9  %15
    0.85 %16
    ];
Node_Demand=[
    0   %1
    30  %2 30
    90  %3 90
    0   %4
    0   %5
    0   %6
    0   %7
    0   %8
    0   %9
    0   %10
    0   %11
    110 %12 110
    40 %13 40
    90  %14 90
    ];

Compressor_Location=[%Node Efficiency Ca,Cb,Cc
    5
    7
    9
    11
    ];
Compressor_Efficiency=[
    0.83
    0.84
    0.83
    0.84
    ];
Compressor_Costs=[
    0 0.2e-3 0
    0 0.2e-3 0
    0 0.2e-3 0
    0 0.2e-3 0
    ];
Compressor_Ratio=[
    1.6 1
    1.8 1
    1000 2
    1100 2
    ];

%Standard temperature
Standard_Temp=492.3;%R
%Standard pressure
Standard_Pressure=14.5038;%psia
%Average gas temperature
Average_Temperature=520;%R
%Gas specific gravity
Gas_Gravity=0.6;%no unit
%Gas compressibility factor
Gas_Comp=0.87;%9987;%No units
%Compressor suction temperature
Compressor_Temp=520;
%Gas specific heat ratio
Gas_Heat=1.3049;
%Gas compressibility factor at compressor inlet
Compressor_iComp=0.9987;



