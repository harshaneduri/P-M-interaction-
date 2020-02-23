%% 

clc;
clear;
%   Input parameters  %
b=300;      % Least dimension of column
D=600;      % Greater dimension of column
fck=65;     % Grade of concrete
fy=415;     % Grade of steel
lyr=3;      % No. of layers of steel
cnt=lyr;
k_count=3;
%%
% Theoretical maximum axial compression %
% p_steel=[0.8,1.2,1.6,2];
% pbyfck=0.1;
p_steel=0.8;
depth_ratio=[0.05,0.1,0.15,0.2];
% plotName1="P-M interaction curve";
% pstar=["p%=";"p%=";"p%=";"p%="];
% delta=["dc/D=";"dc/D=";"dc/D=";"dc/D="];
% plotName2=[pstar,num2str(p_steel')];
% plotName3=[delta,num2str(depth_ratio')];
% if fck==65
%     e1=0.0021;
%     e2=0.0033;
%     n1=1.9;
% elseif fck==70
%     e1=0.0022;
%     e2=0.0031;
%     n1=1.7;
% elseif fck==75
%     e1=0.0023;
%     e2=0.0029;
%     n1=1.6;
% elseif fck==80
%     e1=0.0023;
%     e2=0.0028;
%     n1=1.5;
% elseif fck==85
%     e1=0.0024;
%     e2=0.0027;
%     n1=1.5;
% elseif fck==90
%     e1=0.0024;
%     e2=0.0026;
%     n1=1.4;
% end

Es=2e5;
if (fy>250)
    eps_sy=(0.87*fy/Es)+0.002;
else
    eps_sy=(0.87*fy/Es);
end
for m=1
    dd=depth_ratio(m)*D;
    dis(1,1)=D-dd; dis(2,1)=0.5*D; dis(3,1)=dd;
    depth_ratio(m)
    for l=1
        As(1,1)=p_steel(l)*b*D/300;
        As(2,1)=p_steel(l)*b*D/300;
        As(3,1)=p_steel(l)*b*D/300;
        p_steel(l)
        
        for j=1
            %IS formula
%             if j==1
%                 eps_cp=0.002;
%                 eps_cu=0.0035;
%                 n=2;
%             end
            %67.5 mpa concrete values
%             if j==2
%                 eps_cp=0.0022;
%                 eps_cu=0.0032;
%                 n=1.793;
%             end
           
%                 eps_cp=1.357*fck/10^5 +0.001225;
%                 eps_cu=-3*fck/10^5 + 0.0052;
%                 n=-0.02*fck + 3.1571;
 %77.5 mpa concrete values
%             if j==2
%                 eps_cp=0.0023;
%                 eps_cu=0.0028;
%                 n=1.544;
%             end
%             %87.5 mpa concrete values
%             if j==2
%                 eps_cp=0.0024;
%                 eps_cu=0.0027;
%                 n=1.4374;
%             end
            if j==1
               eps_cp=(2+0.085*((0.8*fck-50)^0.53))/10^3;
                eps_cu=(2.6+35*((90-0.8*fck)/100)^4)/10^3;
                n=1.4+23.4*((90-0.8*fck)/100)^4 ;
            end
%             if j==6
%                 eps_cp=e1;
%                 eps_cu=e2;
%                 n=n1;
% %             end
            if eps_cp>=0.0024
                Puo=(0.453*fck*b*D)+(0.8265*fy-0.453*fck)*(As(1,1)+As(2,1)+As(3,1));
            elseif eps_cp>=0.0023
                Puo=(0.453*fck*b*D)+(0.8163*fy-0.453*fck)*(As(1,1)+As(2,1)+As(3,1));
            elseif eps_cp>=0.0022
                Puo=(0.453*fck*b*D)+(0.8074*fy-0.453*fck)*(As(1,1)+As(2,1)+As(3,1));
            elseif eps_cp>=0.0021
                Puo=(0.453*fck*b*D)+(0.7986*fy-0.453*fck)*(As(1,1)+As(2,1)+As(3,1));
            elseif eps_cp>=0.002
                Puo=(0.453*fck*b*D)+(0.7897*fy-0.453*fck)*(As(1,1)+As(2,1)+As(3,1));
            end
%             Puo=(0.453*fck*b*D)+(0.7897*fy-0.453*fck)*(As(1,1)+As(2,1)+As(3,1));
            i=1;
            increments=40;clean=1;indexing=0;
            for P=0:Puo/increments:((increments-1)/increments)*Puo
                e=1;iter=0;
                xu=0.5*D;
                while (e>0.001)        
%                 for xu=dd:0.1:15*D 
                    St=0;
                    Sc=0;
                    MSt=0;
                    MSc=0;
                    if xu<=D
                        eps_c=eps_cu;
                    else
                        y=D*eps_cp/eps_cu;
                        eps_c=eps_cp*((1+((D-y)/(xu-(D-y)))));
                        eps_cmin=eps_c*(xu-D)/xu;
                    end
                    for count=1:cnt
                        eps_st(count,1)=eps_c*(dis(count,1)-xu)/xu;
                        xi=abs(eps_st(count,1));
                        if xi>eps_sy
                            yi=0.87*fy;
                        elseif xi>=(0.87*0.975*fy/Es)+0.001
                                yi=((xi-((0.87*0.975*fy/Es)+0.001))*(0.87*fy-0.84825*fy)/(eps_sy-((0.87*fy/Es)+0.001)))+0.84825*fy;
                        elseif xi>=(0.87*0.95*fy/Es)+0.0007
                                yi=((xi-((0.87*0.95*fy/Es)+0.0007))*(0.84825*fy-0.8265*fy)/(((0.87*fy/Es)+0.001)-((0.87*0.95*fy/Es)+0.0007)))+0.8265*fy;
                        elseif xi>=(0.87*0.9*fy/Es)+0.0003
                                yi=((xi-((0.87*0.9*fy/Es)+0.0003))*(0.8265*fy-0.783*fy)/(((0.87*0.95*fy/Es)+0.0007)-((0.87*0.9*fy/Es)+0.0003)))+0.783*fy;
                        elseif xi>=(0.87*0.85*fy/Es)+0.0001
                                yi=((xi-((0.87*0.85*fy/Es)+0.0001))*(0.783*fy-0.7395*fy)/(((0.87*0.9*fy/Es)+0.0003)-((0.87*0.85*fy/Es)+0.0001)))+0.7395*fy;
                        elseif xi>=(0.87*0.8*fy/Es)
                                yi=((xi-(0.87*0.8*fy/Es))*(0.7395*fy-0.696*fy)/(((0.87*0.85*fy/Es)+0.0001)-(0.87*0.8*fy/Es)))+0.696*fy;
                        elseif xi>=0
                                yi=(xi*Es);
                        elseif xi<0
                                yi=0;
                        end
                        fst(count,1)=yi;
                    end
                    if (xu<=D)
                        fc1=@(e) fck*(1-(1- e/eps_cp).^n);
                        fc2=@(e) fck*e.^(0);
                        a=0.447*xu*(integral(fc1,0,eps_cp)+integral(fc2,eps_cp,eps_c))/(fck*eps_c)/D;
                        fc1_cen=@(e) fck*(1-(1-e/eps_cp).^n).*e;
                        fc2_cen=@(e) fck*e;
                        k2=1-(integral(fc1_cen,0,eps_cp)+integral(fc2_cen,eps_cp,eps_c))/(eps_c*(integral(fc1,0,eps_cp)+integral(fc2,eps_cp,eps_c)));
%                         a=0.362*xu/D;
                        xdash=k2*xu;
                    else
%                         g=16/((7*xu/D)-3)^2;
                        g=(y/(xu-(D-y)))^n;
                        a=0.447*(1-g*y/((n+1)*D));
                        fc1=@(e) fck*(1-(1- e/eps_cp).^n);
                        fc1_cen=@(e) fck*(1-(1-e/eps_cp).^n).*e;
                        x_cen=(integral(fc1_cen,eps_cmin,eps_cp)/(eps_c*integral(fc1,eps_cmin,eps_cp)))*xu;
                        xfinal=(fck*y*(xu-D+y/2)-((integral(fc1,eps_cmin,eps_cp))/(eps_cp-eps_cmin))*y*(x_cen))/(fck*y-((integral(fc1,eps_cmin,eps_cp))/(eps_cp-eps_cmin))*y);
                        xdash=(0.5*D*D-(g*y/(n+1))*(xu-xfinal))/(D-(g*y/(n+1)));
%                         xdash=(0.5*D*D-(g*y/3)*(D-y+3*y/4)/(D-(g*y/3)));
                    end
                    for counter=1:cnt
                        if eps_st(counter,1)>0
                            St=St+(fst(counter,1)*As(counter,1));
                            MSt=MSt+(fst(counter,1)*As(counter,1))*((dis(counter)-(0.5*D)));
                        else
                            Sc=Sc+(fst(counter,1)*As(counter,1));
                            MSc=MSc+(fst(counter,1)*As(counter,1))*((0.5*D)-(dis(counter)));
                        end
                    end
                    xu_calc=(P+St-Sc)/(a*fck*b*D/xu);
                    e=abs(xu_calc-xu)/xu;
                    if (e<=0.001)
                          break;
                    else
                        xu_new=0.5*(xu+xu_calc);
                    end  
                    xu=xu_new;
                    iter=iter+1;
                    if iter>300
                        break;
                    end
                end
                
                C=a*fck*b*D;
                if iter<300
                    M(i,j)=(C*(0.5*D-xdash)+MSc+MSt);
                    Load(i,j)=P;
                    
                else
                    indexing(clean)=i;
                    clean=clean+1;
                end 
                i=i+1;
%                 disp(i)
            end
            
            M(i,j)=0;
            Load(i,j)=Puo;
            disp(j)
            if any(indexing==0)
            else
                for k=1:length(indexing)
                    M(indexing(k),j)=M(indexing(k)-1,j);
                    Load(indexing(k),j)=Load(indexing(k)-1,j);
                end
            end
%             if indexing(j)<increments+1
%                 M(indexing(j)+1:increments+1,j)=M(indexing(j),j);
%                 Load(indexing(j)+1:increments+1,j)=Load(indexing(j),j);
%             end
        end
%         min_ind=min(indexing);
%         max_ind=max(indexing);     
%         for j=1:6
%             if indexing(j)==max_ind
%             else
%                 M(indexing(j):max_ind-1,j)= M(indexing(j)-1,j);
%                 Load(indexing(j):max_ind-1,j)= Load(indexing(j)-1,j);
%             end
%         end
%         figure
%             if l<=2
%                 subplot(1,2,l);
%             else
%                 subplot(1,2,(l-2));
% %             end
        plot(M(:,1)/(fck*b*D*D),Load(:,1)/(fck*b*D),'LineWidth',2,"Color",'r');
%             hold on 
%             plot(M(:,2)/(fck*b*D*D),Load(:,2)/(fck*b*D),'LineWidth',2,"Color",'b');
%             plot(M(:,3)/(fck*b*D*D),Load(:,3)/(fck*b*D),'LineWidth',2,"Color",'k');
%             plot(M(:,4)/(fck*b*D*D),Load(:,4)/(fck*b*D),'LineWidth',2,"Color",'g');
%             plot(M(:,5)/(fck*b*D*D),Load(:,5)/(fck*b*D),'LineWidth',2,"Color",'m',"Marker",'o',"MarkerSize",3);
% %             plot(M(:,6)/(fck*b*D*D),Load(:,6)/(fck*b*D),'LineWidth',2,"Color",'k');
%         axis([0,inf,0,inf])
%         xlabel('M/(f_{ck}bD^{2})');
%         ylabel('P/(f_{ck}bD)');
%         name=[char(plotName1),strcat(plotName2(l,1),plotName2(l,2)),strcat(plotName3(m,1),plotName3(m,2))];
%         title(name);
%         legend('IS formula','Average values (67.5 MPa)','IRC 112 formula',"Location","southoutside");
%         hold off
        Ld=Load/(fck*b*D);
        Md=M/(fck*b*D*D);
        
        A=table(Md(:,1),Ld(:,1));    
%         B=table(Ld(:,2),Md(:,2));
%         C=table(Ld(:,3),Md(:,3));
% %         k_count=2;
%         writetable(A,"HSCv3.xlsx","sheet",k_count,"Range",'E47');
%         writetable(B,"trail.xlsx","sheet",k_count,"Range",'D47');
%         writetable(C,"trail.xlsx","sheet",k_count,"Range",'G47');
%         writetable(B,"HSC_values.xlsx","sheet",1,"Range",'B3');
%         MISdeviation(m,l)=abs(mean(M(1:end-1,1)./M(1:end-1,3)));
%         M_Avg_deviation(m,l)=abs(mean(M(1:end-1,2)./M(1:end-1,3)));
%         M_A2_deviation(m,l)=abs(mean(M(1:end-1,3)./M(1:end-1,6)));
%         M_A3_deviation(m,l)=abs(mean(M(1:end-1,4)./M(1:end-1,6)));
%         M_IRC_deviation(m,l)=abs(mean(M(1:end-1,5)./M(1:end-1,6)));
    end   
end

% psteel_08=MISdeviation(:,1);
% psteel_12=MISdeviation(:,2);
% psteel_16=MISdeviation(:,3);
% psteel_20=MISdeviation(:,4);
% depth_ratios=depth_ratio';
% M_IS_deviation=table(depth_ratios,psteel_08,psteel_12,psteel_16,psteel_20);
% psteel_08=Mavgdeviation(:,1);
% psteel_12=Mavgdeviation(:,2);
% psteel_16=Mavgdeviation(:,3);
% psteel_20=Mavgdeviation(:,4);
% M_avg_deviation=table(depth_ratios,psteel_08,psteel_12,psteel_16,psteel_20);
% psteel_08=Minterpolatedeviation(:,1);
% psteel_12=Minterpolatedeviation(:,2);
% psteel_16=Minterpolatedeviation(:,3);
% psteel_20=Minterpolatedeviation(:,4);
% M_interpolate_deviation=table(depth_ratios,psteel_08,psteel_12,psteel_16,psteel_20);
% writetable(M_IS_deviation,"meandeviationvalues.xlsx","Sheet",1,"Range",'A1');
% writetable(M_avg_deviation,"meandeviationvalues.xlsx","Sheet",2,"Range",'A1');
% writetable(M_interpolate_deviation,"meandeviationvalues.xlsx","Sheet",3,"Range",'A1');