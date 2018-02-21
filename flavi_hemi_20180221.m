function [  ] = monteTri(howLong, pH, Pdim)
%modified version of I-H Lee's Pop-out model (2011) for modeling flavi hemifusion
%Each monomer is at equilibrium between an extended ('X' and retracted 'R'
%state) - only the extended state can trimerize (go to state 'T')
%monomers form a 30-mer surface patch, where adjacent protomers can
%influence their n+1 neighbor (only their right side neighbors within a dimer).
%Extended subunit increases the activation probability of adjacent subunit.
%Trimers require n+1 and n+2 subunits to also extend.
%Final step when two trimers are formed (two occurances of n+1 and n+2
%extended. Remaining unextended are deleted.
%Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

t = 0;
step = 1; %ms
plotStep = 1; %ms unit time for recording the history

global H;
H = 10^(-pH);

numOfSteps = floor(howLong/step);
timeHistory = 0;

virions = 500;
Etot=virions*(30/(6.022 * 10^23)); %total concentration E (scales with number of virions)
Es=zeros(virions,30);
Edist=zeros(1,5);
EdistHistory=[virions*30 0 0 0 0];

%state vector for flavi virion - 30 subunits, each with states
%0: retracted
%1: extended
%2: trimerized
%3: zipped (hemifused)
%4: dead subunit (incompetent for trimerization)
%structure do be defined below

Hhistory = [0];

detMonteIndicator = 0; % variable to keep balance coming from deterministic monte carlo coupling

P = zeros(5,2);
% Probability matrix for MC simulation
dice = 0;
% a variable to store generated random #

k_act = 1 * H; %rate constant for extension (mM-1 ms-1)
k_ret = 3.1623 * 10^-6; %rate constant for retraction (ms-1)
K_DiMono = k_act/k_ret; %
k_dead = 0;
k_tri = .01; %rate constant for trimerization
k_tri_dis = 1 * 10^-9; %rate constant for trimer dissociation
k_zip = .5;
dimer_dissFactor = Pdim; %amount to increase extension probability by neighbors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%initialization of the system using E dimer-monomer equilibirium%%%%
for j=1:virions
    for k=1:30
        dice = rand;
        if dice < (K_DiMono)/(K_DiMono+1)
            Es(j,k) = 1;
        end
        if dice < k_dead
            Es(j,k) = 4;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:numOfSteps
    
    if mod(i*step,plotStep) == 0
        timeHistory = [timeHistory;i*step];
    end
    
    if rand() < 0.5
        detMonteIndicator = 1;
    else
        detMonteIndicator = 0;
    end
 
    %%%%%Setting up probabilities for reactions - where the different rates are set%%%%%%%%
    P(1,1) = 1-exp(-k_act * step);
    P(1,2) = 1-exp(-k_ret * step);
    P(1,3) = 1-exp(-k_dead * step);
    P(2,1) = 1-exp(-k_tri * step);
    P(2,2) = 1-exp(-k_tri_dis * step);
    P(3,1) = 1-exp(-k_zip * step);
    P(3,2) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%setting up reactions in the structure
    for j=1:1:virions
        k = floor(rand()*30) +1 ;
        for k_index =1:1:30
            k = k + k_index;
            if k > 30
                k = k - 30;
            end
            dice = rand();
            switch Es(j,k)
                %set a certain proportion dead
                case 0
                    %dimer dissociation step
                     if k == 1
                        if Es(j,30) == 1 || Es(j,2) == 1
                            if dice < ( 1-exp(-k_act * step * Pdim) )
                                Es(j,k) = 1;
                            end
                        else
                            if dice < P(1,1)
                                Es(j,k) = 1;
                            end
                        end
                     elseif k == 30
                        if (Es(j,1) == 1) || (Es(j,29) == 1)
                            if dice < ( 1-exp(-k_act * step * Pdim) )
                                Es(j,k) = 1;
                            end
                        else
                            if dice < P(1,1)
                                Es(j,k) = 1;
                            end
                        end
                     elseif Es(j,k+1) == 1
                        if (mod(k,2) ~= 0)
                            if dice < ( 1-exp(-k_act * step * Pdim) )
                                Es(j,k) = 1;
                            end
                        else
                            if dice < P(1,1)
                                Es(j,k) = 1;
                            end
                        end   
                     else
                        if Es(j,k+1) == 1 && Es(k,k+2) == 1 && (mod(k,2) ~= 0)
                            if dice < (( 1-exp(-k_act * step * Pdim))*1.5)
                                Es(j,k) = 1;
                            end
                        else
                            if dice < P(1,1)
                                Es(j,k) = 1;
                            end
                        end
                     end
                case 1
                    %trimerization only when two adjacent subunits are
                    %extended
                    if k <= 28
                        if (k ~= 6) && (k ~= 7) && (k ~= 17) && (k ~= 18) && (k ~= 26) && (k ~= 27) && (Es(j,k+1) == 1) && (Es(j,k+2) == 1)
                            if dice < ( P(2,1) )
                                Es(j,k+1) = 2;
                                Es(j,k+2) = 2;
                                Es(j,k) = 2;
                            end
                        else
                            if dice < ( P(1,2) )
                                Es(j,k) = 0;
                            end
                        end
                    elseif k <= 12
                        if (((k == 1) || (k == 3) || (k == 5) || (k == 7)) && (Es(j,k+8) == 1) && (Es(j,k+10) == 1))
                            if dice <  ( P(2,1) )
                                Es(j,k+8) = 2;
                                Es(j,k+10) = 2;
                                Es(j,k) = 2;
                            end
                        else
                            if dice < ( P(1,2) )
                                Es(j,k) = 0;
                            end
                        end
                    elseif ((k == 10) || (k == 12) || (k == 14) || (k == 16))
                        if (Es(j,k+9) == 1) && (Es(j,k+11) == 1)
                            if dice <  ( P(2,1) )
                                Es(j,k+9) = 2;
                                Es(j,k+11) = 2;
                                Es(j,k) = 2;
                            end
                        else
                            if dice < ( P(1,2) )
                                Es(j,k) = 0;
                            end
                        end
                    elseif (k == 20) || (k == 22) || (k == 24)
                        if (Es(j,k+2) == 1) && (Es(j,k+7) == 1)
                            if dice <  ( P(2,1) )
                                Es(j,k+2) = 2;
                                Es(j,k+7) = 2;
                                Es(j,k) = 2;
                            end
                        else
                            if dice < ( P(1,2) )
                                Es(j,k) = 0;
                            end
                        end      
                    else
                        if (((k == 8) || (k == 10) || (k == 12) || (k == 14) || (k == 16)) && (Es(j,k+2) == 1) && (Es(j,k+11) == 1))
                            if dice < ( P(2,1) )
                                Es(j,k) = 2;
                            end
                        else
                            if dice < ( P(1,2) )
                                Es(j,k) = 0;
                            end
                        end
                        %else
                        %    if dice < ( P(1,2) )
                        %        Es(j,k) = 0;
                        %    end
                    end
                case 2
                    %zipping (from trimer to zipped, or trimer back to
                    %active monomers)
                    if (k == 1) || (k == 8) || (k == 11) || (k == 19) || (k == 22)
                        if (Es(j,k+1) == 2 && Es(j,k+2) == 2 && Es(j,k+3) == 2 && Es(j,k+4) == 2 && Es(j,k+5) == 2)
                            if dice < ( P(3,1) )
                                %Es(j,k) = 3;
                                %Es(j,k+1) = 3;
                                %Es(j,k+2) = 3;
                                %Es(j,k+3) = 3;
                                %Es(j,k+4) = 3;
                                %Es(j,k+5) = 3;
                                Es(j,:)=3;
                            end
                        else
                            if dice < ( P(2,2) )
                                Es(j,k) = 1;
                            end
                        end
                    elseif (k == 1) || (k == 4)
                        if (Es(j,k+1) == 2 && Es(j,k+2) == 2 && Es(j,k+7) == 2 && Es(j,k+8) == 2 && Es(j,k+9) == 2)
                            if dice < ( P(3,1) )
                                %Es(j,k) = 3;
                                %Es(j,k+1) = 3;
                                %Es(j,k+2) = 3;
                                %Es(j,k+7) = 3;
                                %Es(j,k+8) = 3;
                                %Es(j,k+9) = 3;
                                Es(j,:)=3;
                            end
                        else
                            if dice < ( P(2,2) )
                                Es(j,k) = 1;
                            end
                        end
                    elseif (k == 11) || (k == 14)
                        if (Es(j,k+1) == 2 && Es(j,k+2) == 2 && Es(j,k+8) == 2 && Es(j,k+9) == 2 && Es(j,k+10) == 2)
                            if dice < ( P(3,1) )
                                %Es(j,k) = 3;
                                %Es(j,k+1) = 3;
                                %Es(j,k+2) = 3;
                                %Es(j,k+8) = 3;
                                %Es(j,k+9) = 3;
                                %Es(j,k+10) = 3;
                                Es(j,:)=3;
                            end
                        else
                            if dice < ( P(2,2) )
                                Es(j,k) = 1;
                            end
                        end    
                    elseif (k == 1) || (k == 10)
                        if (Es(j,k+2) == 2 && Es(j,k+4) == 2 && Es(j,k+8) == 2 && Es(j,k+10) == 2 && Es(j,k+12) == 2)
                            if dice < ( P(3,1) )
                                %Es(j,k) = 3;
                                %Es(j,k+2) = 3;
                                %Es(j,k+4) = 3;
                                %Es(j,k+8) = 3;
                                %Es(j,k+10) = 3;
                                %Es(j,k+12) = 3;
                                Es(j,:)=3;
                            end
                        else
                            if dice < ( P(2,2) )
                                Es(j,k) = 1;
                            end
                        end
                    elseif (k == 3) || (k == 12)
                        if (Es(j,k+2) == 2 && Es(j,k+4) == 2 && Es(j,k+10) == 2 && Es(j,k+12) == 2 && Es(j,k+14) == 2)
                            if dice < ( P(3,1) )
                                %Es(j,k) = 3;
                                %Es(j,k+2) = 3;
                                %Es(j,k+4) = 3;
                                %Es(j,k+10) = 3;
                                %Es(j,k+12) = 3;
                                %Es(j,k+14) = 3;
                                Es(j,:)=3;
                            end
                        else
                            if dice < ( P(2,2) )
                                Es(j,k) = 1;
                            end
                        end
                    else
                        if (k == 22) && (Es(j,k+2) == 2 && Es(j,k+4) == 2 && Es(j,k+6) == 2 && Es(j,k+7) == 2 && Es(j,k+8) == 2)
                             if dice < ( P(3,1) )
                                %Es(j,k) = 3;
                                %Es(j,k+2) = 3;
                                %Es(j,k+4) = 3;
                                %Es(j,k+6) = 3;
                                %Es(j,k+7) = 3;
                                %Es(j,k+8) = 3;
                                Es(j,:)=3;
                             end
                        else
                             if dice < ( P(2,2) )
                                Es(j,k) = 1;
                             end
                        end
                        %else
                        %    if dice < P(2,2)
                        %        Es(j,k) = 1;
                        %    end
                    end
                %un-zipping
                %case 3
                %    if dice < (P(3,2))
                %        Es(j,k) = 2;
                %    end
                    
            end
        end
    end

    %this below sums the number of each species across all 30 Es
    Edist= hist(Es(:,1),[0 1 2 3 4])+hist(Es(:,2),[0 1 2 3 4])+hist(Es(:,3),[0 1 2 3 4])+hist(Es(:,4),[0 1 2 3 4])+hist(Es(:,5),[0 1 2 3 4])+hist(Es(:,6),[0 1 2 3 4])+hist(Es(:,7),[0 1 2 3 4])+hist(Es(:,8),[0 1 2 3 4])+hist(Es(:,9),[0 1 2 3 4])+hist(Es(:,10),[0 1 2 3 4])+hist(Es(:,11),[0 1 2 3 4])+hist(Es(:,12),[0 1 2 3 4])+hist(Es(:,13),[0 1 2 3 4])+hist(Es(:,14),[0 1 2 3 4])+hist(Es(:,15),[0 1 2 3 4])+hist(Es(:,16),[0 1 2 3 4])+hist(Es(:,17),[0 1 2 3 4])+hist(Es(:,18),[0 1 2 3 4])+hist(Es(:,19),[0 1 2 3 4])+hist(Es(:,20),[0 1 2 3 4])+hist(Es(:,21),[0 1 2 3 4])+hist(Es(:,22),[0 1 2 3 4])+hist(Es(:,23),[0 1 2 3 4])+hist(Es(:,24),[0 1 2 3 4])+hist(Es(:,25),[0 1 2 3 4])+hist(Es(:,26),[0 1 2 3 4])+hist(Es(:,27),[0 1 2 3 4])+hist(Es(:,28),[0 1 2 3 4])+hist(Es(:,29),[0 1 2 3 4])+hist(Es(:,30),[0 1 2 3 4]);
    %below generates data binned per virion
    %Edist= hist(Es(1,:),[0 1 2 3])+hist(Es(2,:),[0 1 2 3])+hist(Es(3,:),[0 1 2 3])+hist(Es(4,:),[0 1 2 3])+hist(Es(5,:),[0 1 2 3])+hist(Es(6,:),[0 1 2 3])+hist(Es(7,:),[0 1 2 3])+hist(Es(8,:),[0 1 2 3])+hist(Es(9,:),[0 1 2 3])+hist(Es(10,:),[0 1 2 3])
    %Ehist= hist(Es(1,:),[3])
    %Enum = nnz(Es(1:virions,:)==3)
    %hist(
        
    if mod(i*step,plotStep) == 0
        EdistHistory = [EdistHistory;Edist];
        %HemiHist = [EdistHistory;Hemi];
        %Hemi = unique(Es(1:1:virions,:)==1);
    end
    
    if H < 0
        H = 0;
    
    if mod(i,100) == 0
        i;
    end
end
%csvwrite(['pH ',num2str(pH),'k_act ',num2str(k_act),'k_dis ',num2str(k_dis),'dimer_diss_Factor ',num2str(Pdim),'k_tri ',num2str(k_tri),'k_tri_dis ',num2str(k_tri_dis),'k_zip ',num2str(k_zip),'Param.txt'],EdistHistory);
csvwrite('file3.csv',EdistHistory)

%figure(1); clf;
%figure(2); clf;
plot(EdistHistory);
%plot(hist(EdistHistory));
%plot(hist(Hemi));
%plot(Enum);
%plot(hist(Hemi,timehistory))
%figure(3); clf;

toc;
    
end