clear all
close all
clc
pause on
maxDistance=60;
samples=1;
numSolved=0;
numRelevant=0;
NUMCONNECTIONS=1;
numTopologies=1;
var_vector=2;%[0 2 4];
for m=1:length(var_vector)
    a_s = [sqrt(var_vector(m))]; %Used in the beamforming weight design.
    for topologyID=0:numTopologies-1
        %Secondary network topology
        xVector=[4.25 61.551762 60.516459 13.727010 20.250175 30.520528 63.421532 90.844299 91.446 75.342];
        yVector=[45.884092 20.654070 85.410185 75.017504 69.373124 20.183925 50.813904 23.001942 58.42 93.245];
        numNodes=length(xVector);
        %Position of primary receiver
        xPrimary=33.616059;
        yPrimary=51.549543;
        xVectorTot=[xVector xPrimary];
        yVectorTot=[yVector yPrimary];
        
        EP_SINR_matrix=zeros(numNodes,numNodes);
        EP_SINR_matrix_noOpt=zeros(numNodes,numNodes);
        ES_SINR_matrix=zeros(numNodes,numNodes);
        ES_SINR_matrix_noOpt=zeros(numNodes,numNodes);
        positionPlotting(xVectorTot, yVectorTot, 11, 3,0,100,0,100)
        %loadConnectionList;
        %loadConnectionLengths;
        relevantIRFvaluesSum=0;
        relevantIRFvaluesNumber=0;
        interferenceReductionFactor=ones(length(xVector),length(xVector),length(xVector));
        for i=1:length(xVector)
            for j=1:length(xVector)
                fprintf(1,'Considering transmitter N%d to receiver N%d\n',i-1,j-1);
                %pause
                d=distance(xVector,yVector,i,j);
                if((i~=j)&&(d<=maxDistance))&&(i==7)&&(j==9)
                    
                    [thetap Kp affectedNodesIndexes]=evaluate_DOA_Vector(xVectorTot,yVectorTot,i,j,maxDistance)
                    K = 1; %1 intended target
                    %Kp=1;%Kp devices to be avoided
                    Mt=8; %Number of transmit Antennas at the Transmitter
                    Mr=1; %1 Receive Antenna per user
                    
                    %theta=acosd((xVector(j)/d));
                    
                    if(xVector(i)>xVector(j))
                        theta=180+atand((yVector(i)-yVector(j))/(xVector(i)-xVector(j)));
                    elseif (yVector(i)>yVector(j))
                        theta = 360+atand((yVector(i)-yVector(j))/(xVector(i)-xVector(j)));
                    else
                        theta=atand((yVector(j)-yVector(i))/(xVector(j)-xVector(i)));
                    end
%                     if(yVector(i)<yVector(j))
%                         theta=360-theta;
%                     end
                    %theta=[35 60];%AOD for the 2 Secondary Users (measured)
                    %thetap=;%AOD for the Primary User
                    intfr = 0.05;%Maximum Allowable Interference to Primary per Secondary User
                    
                    warning off
                    test=[0 1 2 3 4 5 6 7];%Antenna Element Number
                    
                    %cx = 4;
                    %            for check =1:1:length(a_s)
                    
                    consr2=ones(samples,length(thetap));
                    consr2_improvement_ratio=ones(samples,length(thetap));
                    consr2_noOpt=ones(samples,length(thetap));
                    for cx = 1:1:samples
                        numRelevant=numRelevant+1;
                        [R]=covam(theta,a_s,Mt);% Covariance Secondary User
                        for q=1:length(thetap)
                            [Rp(:,:,q)] =covam(thetap(q),a_s,Mt);%Covariance for Primary and secondary affected users
                        end
                        %***Solving Minimum Transmission Power Problem with interference constraint***
                        cvx_solver sedumi
                        cvx_begin SDP
                        cvx_precision high
                        cvx_quiet(true)
                        variable g(Mt,Mt) complex
                        minimize(abs(trace(g)))
                        subject to
                        trace(R*g)-10*0.01==semidefinite(1);
                        for q=1:length(thetap)
                            intfr-trace(Rp(:,:,q)*g)==hermitian_semidefinite(1);
                        end
                        
                        
                        g==g'
                        g>=0;
                        cvx_end
                        if strcmp(cvx_status,'Solved')
                            numSolved=numSolved+1;
                            %********Saving Results*******
                            power2(cx)=cvx_optval; %Optimal minimised power allocation
                            
                            snr1(cx)=(trace(R*g))/0.01;                        %snr at the measured location of Sec user 1
                            c=0.01;
                            %                     c_old=c;
                            %                     c_min=c;
                            %                     c_max=c;
                            c_step=1.1;
                            snr1_noOpt(cx)=trace(R*(c*ones(Mt,Mt)))/0.01;   %snr at the measured location of Sec user 1
                            oldSign=sign(real(snr1(cx)-snr1_noOpt(cx)));
                            oldValue=real(snr1(cx)-snr1_noOpt(cx));
                            while abs(real(snr1(cx)-snr1_noOpt(cx)))>2
                                if(real(snr1(cx)-snr1_noOpt(cx))>0)
                                    c=c_step*c;
                                else
                                    c=c/c_step;
                                end
                                snr1_noOpt(cx)=trace(R*(c*ones(Mt,Mt)))/0.01;
                                %c;
                                c_step;
                                newValue=abs(real(snr1(cx)-snr1_noOpt(cx)));
                                newSign=sign(real(snr1(cx)-snr1_noOpt(cx)));
                                fprintf(1,'newValue: %f; oldValue: %f: newSign: %d; oldSign: %d; difference: %f\n',newValue, oldValue, newSign, oldSign, real(snr1(cx)-snr1_noOpt(cx)));
                                if newSign~=oldSign
                                    if newValue<=oldValue
                                        break
                                    else
                                        oldSign=newSign;
                                    end
                                else
                                    if abs(newValue-oldValue)<1e-6
                                        c_step_old=c_step
                                        c_step=0.9+0.1*c_step+(c_step -(0.9+0.1*c_step))*rand(1,1)
                                        if c_step==c_step_old
                                            break
                                        end
                                    end
                                end
                                oldValue=newValue;
                            end
                            for q=1:length(thetap)
                                consr2(cx,q)=trace(Rp(:,:,q)*g);  %Interference at the primary user.
                                consr2_noOpt(cx,q)=trace(Rp(:,:,q)*(c*ones(Mt,Mt)));  %Interference at the primary user without optimization.
                                consr2_improvement_ratio(cx,q)=consr2_noOpt(cx,q)/consr2(cx,q);
                            end
                            r_loc1 = theta + normrnd(0,a_s);                        %real location of Sec user 1, normally distributed
                            [R1]=covam(r_loc1,0,Mt);
                            SINR1(cx) = trace(R1*g)/0.01;                         %SNR at the actual position of secondary user 1
                            SINR1_noOpt(cx) = trace(R1*(c*ones(Mt,Mt)))/0.01;      %SNR at the actual position of secondary user 1 without optimization
                            cx;
                        else
                            snr1(cx)=10;
                            snr1_noOpt(cx)=10;
                            SINR1(cx)=10;
                            SINR1_noOpt(cx)=10;
                            
                        end
                    end
                    
                    
                    if(samples>1)
                        consr2_improvement_ratio_average=mean(consr2_improvement_ratio);
                    else
                        consr2_improvement_ratio_average=consr2_improvement_ratio;
                    end
                    
                    if length(consr2_improvement_ratio_average)~=length(affectedNodesIndexes)
                        fprintf(1,'Attention\n');
                        fprintf(1,'Something is wrong\n');
                    end
                    interferenceReductionFactor(i,j,affectedNodesIndexes)=real(consr2_improvement_ratio_average);
                    relevantIRFvaluesNumber=relevantIRFvaluesNumber+1;
                    relevantIRFvaluesSum=relevantIRFvaluesSum+mean(real(consr2_improvement_ratio_average));
                    
                    EP_SINR_matrix(i,j)=mean(SINR1);
                    ES_SINR_matrix(i,j)=mean(snr1);
                    
                    EP_SINR_matrix_noOpt(i,j)=mean(SINR1_noOpt);
                    ES_SINR_matrix_noOpt(i,j)=mean(snr1_noOpt);
                    
                    
                    
                end
            end
        end
        
        IRF_average(m)=relevantIRFvaluesSum/relevantIRFvaluesNumber;
        
        EP_SINR_average(m)=mean(mean(EP_SINR_matrix));
        EP_SINR_average_noOpt(m)=mean(mean(EP_SINR_matrix_noOpt));
        
        ES_SINR_average(m)=mean(mean(ES_SINR_matrix));
        ES_SINR_average_noOpt(m)=mean(mean(ES_SINR_matrix_noOpt));
        
        
        consr2_average(m)=mean(mean(consr2));%We take average over samples AND
        %affected terminals: no other way to do it, as number of affected
        %terminals varies for each (i,j) pair
        consr2_average_noOpt(m)=mean(mean(consr2_noOpt)); %We take average over samples AND
        %affected terminals: no other way to do it, as number of affected
        %terminals varies for each (i,j) pair
        %%Output to file to be used as an input for OMNeT++ simulations
        for i=1:length(xVector)
            fileNameStr=sprintf('BF_MUI_Reduction_Factor_error_%d_topology_%d_terminal_%d.txt',var_vector(m),topologyID,i-1);
            omnet_file_ID=fopen(fileNameStr,'a');
            for j=1:length(xVector)
                for k=1:length(xVector)
                    fprintf(omnet_file_ID, '%f ', interferenceReductionFactor(j,k,i));
                end
                fprintf(omnet_file_ID, '\n', interferenceReductionFactor(j,k,i));
            end
            fclose(omnet_file_ID);
        end
        
    end
    fprintf(1, 'Percentage of solved optimizations: %f\n', numSolved/numRelevant);
    
    
    
end


figure()
plot(var_vector,EP_SINR_average);
hold on
plot(var_vector,EP_SINR_average_noOpt,'r');
xlabel('Variance of spread')
ylabel('SNR at real position of intended receiver')
legend('Optimized','Flat')

figure()
plot(var_vector,ES_SINR_average);
hold on
plot(var_vector,ES_SINR_average_noOpt,'r');
xlabel('Variance of spread')
ylabel('SNR at estimated position of intended receiver')
legend('Optimized','Flat')

figure()
plot(var_vector,IRF_average,'LineWidth',2);
xlabel('Variance of spread')
ylabel('Interference Reduction Factor')

% figure()
% plot(var_vector,consr2_average(1,:,1))
% hold on
% plot(var_vector,consr2_average(1,:,2),'r')
% plot(var_vector,consr2_average(1,:,3),'g')
% plot(var_vector,consr2_average(1,:,4),'k')
% xlabel('Variance of spread')
% ylabel('Interference at the position of secondary users and primary system')
%
%
% figure()
% plot(var_vector,consr2_average_noOpt(1,:,1))
% hold on
% plot(var_vector,consr2_average_noOpt(1,:,2),'r')
% plot(var_vector,consr2_average_noOpt(1,:,3),'g')
% plot(var_vector,consr2_average_noOpt(1,:,4),'k')
% xlabel('Variance of spread')
% ylabel('Interference at the position of secondary users and primary system without optimization')