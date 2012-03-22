clear all
close all
pause on
clc
clear all
maxDistance=30;
NUMCONNECTIONS=1;
nn = 100;
DIM_MAX = 200;
choice = 2;
%Secondary network topology
% xVector=[4.25 61.551762 60.516459 13.727010 20.250175 30.520528 63.421532 90.844299 91.446 75.342];
% yVector=[45.884092 20.654070 85.410185 75.017504 69.373124 20.183925 50.813904 23.001942 58.42 93.245];
%Position of primary receiver
% xPrimary=33.616059;
% yPrimary=51.549543;
[dumy,connectionList,xVector,yVector,source,destination] = dijkstras_table(nn,DIM_MAX,choice);
hopLengthList = length(connectionList);
xVectorTot=[xVector];
yVectorTot=[yVector];
positionPlotting(xVectorTot, yVectorTot, 11, 3,0,100,0,100,maxDistance)
% loadConnectionList;
% loadConnectionLengths;
%for m=1:min(NUMCONNECTIONS,length(hopLengthList))
    m = 1;
    for n=1:hopLengthList(m)-1
        [thetap Kp]=evaluate_DOA_Vector(xVectorTot,yVectorTot,connectionList(m,n)+1,connectionList(m,n+1)+1,maxDistance)
        K = 1; %1 intended target
        %Kp=1;%Kp devices to be avoided
        Mt=8; %Number of transmit Antennas at the Transmitter
        Mr=1; %1 Receive Antenna per user
        d=distance(xVector,yVector,connectionList(m,n),connectionList(m,n+1));
%         theta=acosd((xVector(connectionList(i,j+1)+1)/d));
%         if(yVector(connectionList(i,j+1)+1)<yVector(connectionList(i,j)+1))
%             theta=360-theta;
%         end
i=connectionList(m,n);
j=connectionList(m,n+1);
if(xVector(i)>xVector(j))
                        theta=180+atand((yVector(i)-yVector(j))/(xVector(i)-xVector(j)));
                    elseif (yVector(i)>yVector(j))
                        theta = 360+atand((yVector(i)-yVector(j))/(xVector(i)-xVector(j)));
                    else
                        theta=atand((yVector(j)-yVector(i))/(xVector(j)-xVector(i)));
                    end
        %theta=[35 60];%AOD for the 2 Secondary Users (measured)
        %thetap=;%AOD for the Primary User
        intfr = 0.05;%Maximum Allowable Interference to Primary per Secondary User
        samples=2;
        warning off
        test=[0 1 2 3 4 5 6 7]%Antenna Element Number
        var_pos=2:2:6;
        a_s = [sqrt(var_pos)]; %Used in the beamforming weight design.
        %cx = 4;
        for check =1:1:length(a_s)
            for cx = 1:1:samples
                [R(:,:,1)]=covam(theta,a_s(check),Mt);% Covariance Secondary Users
                for q=1:length(thetap)
                    [Rp(:,:,q)] =covam(thetap(q),a_s(check),Mt);%Covariance for Primary
                end
                %***Solving Minimum Transmission Power Problem with interference constraint***
                cvx_solver sedumi
                cvx_begin SDP
                cvx_precision high
                cvx_quiet(true)
                variable g(Mt,Mt,K) complex
                minimize(abs(trace(g(:,:,1))))
                subject to
                for(k=1:1:K)
                    for(l=1:1:K)
                        if(l==k)
                            trace(R(:,:,k)*g(:,:,l))-10*0.01==semidefinite(1);
                            for q=1:length(thetap)
                                intfr-trace(Rp(:,:,q)*g(:,:,1))==hermitian_semidefinite(1);
                            end
                        end
                        if(l~=k)
                            trace(R(:,:,k)*g(:,:,l))==0;
                        end
                    end
                end
                g==g'
                g(:,:,1)>=0;
                cvx_end
                %********Saving Results*******
                power2(cx,check)=cvx_optval; %Optimal minimised power allocation
                
                snr1(cx,check)=(trace(R(:,:,1)*g))/0.01;                        %snr at the measured location of Sec user 1
                c=0.01;
                c_old=c;
                c_min=c;
                c_max=c;
                c_step=1.1;
                snr1_simple(cx,check)=trace(R(:,:,1)*(c*ones(Mt,Mt)))/0.01;   %snr at the measured location of Sec user 1
                oldsign=sign(real(snr1(cx,check)-snr1_simple(cx,check)));
                while abs(real(snr1(cx,check)-snr1_simple(cx,check)))>2
                    if(real(snr1(cx,check)-snr1_simple(cx,check))>0)
                        c=c_step*c
                        %                         c_old
                        %                         if c_old>c
                        %                             c_step=1+0.1*c_step-0.1;
                        %                         end
                    else
                        c=c/c_step
                        %                         c_old
                        %                         if c_old<c
                        %                             c_step=1+0.1*c_step-0.1;
                        %                         end
                    end
                    %                     if c<c_min
                    %                        c_min=c;
                    %                     end
                    %                     if c>c_max
                    %                         c_max=c
                    %                     end
                    %c_old=c;
                    snr1_simple(cx,check)=trace(R(:,:,1)*(c*ones(Mt,Mt)))/0.01;
                    %c;
                    c_step
                    newsign=sign(real(snr1(cx,check)-snr1_simple(cx,check)));
                    real(snr1(cx,check)-snr1_simple(cx,check))
                    if newsign~=oldsign
                        %c_step=0.9+0.1*c_step+(c_step -(0.9+0.1*c_step))*rand(1,1)
                        break
                    else
                        oldsign=newsign;
                    end
                end
                for q=1:length(thetap)
                    consr2(cx,check,q)=trace(Rp(:,:,q)*g);  %Interference at the primary user.
                    consr2_simple(cx,check,q)=trace(Rp(:,:,q)*(c*ones(Mt,Mt)));  %Interference at the primary user without optimization.
                    consr2_improvement_ratio(cx,check,q)=consr2_simple(cx,check,q)/consr2(cx,check,q);
                end
                r_loc1 = theta;% + normrnd(0,sqrt(7));                          %real location of Sec user 1, normally distributed
                [R1]=covam(r_loc1,0,Mt);
                SINR1(cx,check) = trace(R1*g)/0.01;                             %SNR at the actual position of secondary user 1
                SINR1_simple(cx,check) = trace(R1*(c*ones(8,8)))/0.01;      %SNR at the actual position of secondary user 1
                
                cx
                check
            end
        end
        SINR_average=mean(SINR1);
        SINR_average_noOpt=mean(SINR1_simple);
%         figure()
%         plot(var_pos,SINR_average);
%         hold on
%         plot(var_pos,SINR_average_noOpt,'r');
%         xlabel('Variance of spread')
%         ylabel('SNR at real position of intended receiver')
        snr1_average=mean(snr1);
        snr1_average_noOpt=mean(snr1_simple);
%         figure()
%         plot(var_pos,snr1_average);
%         hold on
%         plot(var_pos,snr1_average_noOpt,'r');
%         xlabel('Variance of spread')
%         ylabel('SNR at estimated position of intended receiver')
        
        consr2_average=mean(consr2);
        consr2_average_noOpt=mean(consr2_simple);
%         figure()
%         plot(var_pos,consr2_average(1,:,1))
%         hold on
%         plot(var_pos,consr2_average(1,:,2),'r')
%         plot(var_pos,consr2_average(1,:,3),'g')
%         plot(var_pos,consr2_average(1,:,4),'k')
%         xlabel('Variance of spread')
%         ylabel('Interference at the position of secondary users and primary system')
        
%         figure()
%         plot(var_pos,consr2_average_noOpt(1,:,1))
%         hold on
%         plot(var_pos,consr2_average_noOpt(1,:,2),'r')
%         plot(var_pos,consr2_average_noOpt(1,:,3),'g')
%         plot(var_pos,consr2_average_noOpt(1,:,4),'k')
%         xlabel('Variance of spread')
%         ylabel('Interference at the position of secondary users and primary system without optimization')
        
        
        consr2_improvement_ratio_average=mean(consr2_improvement_ratio);
        consr2_improvement_ratio_average_over_nodes=mean(consr2_improvement_ratio,3);
        interferenceReductionFactor(m,n,:)=consr2_improvement_ratio_average_over_nodes(1,:,1);
    end
%end

