%%% Design of frame
clc
%% Input for the support at the base
s1 = input('Enter 1 or 2 for support hinged or fixed respectively at left side of the frame: ');
s2 = input('Enter 1 or 2 for support roller or hinged respectively at left side of the frame: ');

%% input for beam 
lb = input('Enter the span of the beam: ');

npb = input('Enter the number of the point load acting: ');
mab = 0;
fb = 0;
mbb=0;
for i =1:npb
    pbx(i) = input('Enter the coodinate where point load is acting ');
    pbm(i) = input('Enter the point load magnitute on beam: ');
    mab = mab + (pbx(i)*pbm(i));
    mbb = mbb + ((lb-pbx(i))*pbm(i));
    fb = fb + pbm(i);
end

nub = input('Enter the number of UDL acting on the beam: ');
uc=[0,0,0];
for i = 1:nub
    ubx(i,1) = input('Enter the left end of the UDL coodinate where point load is acting: ');
    ubx(i,2) = input('Enter the right end of the UDL coodinate where point load is acting: ');
    ubm(i) = input('Enter the UDL acting on the beam: ');
    uc = [uc;ubx(i,1),ubx(i,2),ubm(i)];
    mab = mab + (ubm(i)*(ubx(i,2)-ubx(i,1))*0.5*(ubx(i,2)+ubx(i,1)));    
    mbb = mbb + (ubm(i)*(ubx(i,2)-ubx(i,1))*0.5*(lb-(ubx(i,2)+ubx(i,1))));
    fb = fb + (ubm(i)*(ubx(i,2)-ubx(i,1)));
end
uc=uc(2:end,:);
%% Input for left column
disp('Enter the detail of left end of the coulumn:-');
hl = input('Enter the height of the left end of coulumn: ');    

mal = 0;
fl = 0 ;
mbl=0;
malr=0;
npl = input('Enter the number of point load in left coulumn: ');
for i = 1:npl
    plm(i) = input('Enter the magnitude of point load acting: ');
    ply(i) = input('Enter the coordinate where point load is acting: ');
    mal = mal + (plm(i)*ply(i));
    mbl = mbl + (plm(i)*(hl-ply(i)));
    malr = malr  + (plm(i)*(ply(i)-hl+hr));
    fl = fl + plm(i);
end

nul = input('Enter the number of UDL load in left coulumn: ');
for i = 1:nul
    ulm(i) = input('Enter the magnitude of UDL load acting: ');
    uly(i,1) = input('Enter the bottom coordinate where UDL load is started acting: ');
    uly(i,2) = input('Enter the upper coordinate where UDL load is ended acting: ');    
    mal = mal + (ulm(i)*(uly(i,2)-uly(i,1))*0.5*(uly(i,2)+uly(i,1)));
    mbl = mbl + (ulm(i)*(uly(i,2)-uly(i,1))*0.5*(hl-(uly(i,2)+uly(i,1))));
    malr = malr + (ulm(i)*(uly(i,2)-uly(i,1))*(hr-hl+(0.5*(uly(i,2)+uly(i,1)))));
    fl = fl + (ulm(i)*(uly(i,2)-uly(i,1)));
end

%% Input for right column
disp('Enter the detail of right end of the coulumn:-');
hr = input('Enter the height of the right end of coulumn: ');    

 mar = 0;
 fr = 0 ;
 mbr = 0;
 marl=0;
 npr = input('Enter the number of point load in right coulumn: ');
 for i = 1:npr
     prm(i) = input('Enter the magnitude of point load acting: ');
     pry(i) = input('Enter the coordinate where point load is acting: ');
     mar = mar + (prm(i)*pry(i));
     mbr = mbr + (prm(i)*(hr-pry(i)));
     malr = malr  + (plm(i)*(pry(i)+hl-hr));
     fr = fr + prm(i);
 end
 
 nur = input('Enter the number of UDL load in left coulumn: ');
 for i = 1:nur
     urm(i) = input('Enter the magnitude of UDL load acting: ');
     ury(i,1) = input('Enter the bottom coordinate where UDL load is started acting: ');
     ury(i,2) = input('Enter the upper coordinate where UDL load is ended acting: ');    
     mar = mar + (urm(i)*(ury(i,2)-ury(i,1))*0.5*(ury(i,2)+ury(i,1)));
     mbr = mbr + (urm(i)*(ury(i,2)-ury(i,1))*0.5*(hr-(ury(i,2)+ury(i,1))));
     marl = marl + (urm(i)*(ury(i,2)-ury(i,1))*(-hr+hl+(0.5*(ury(i,2)+ury(i,1)))));
    fr = fr + (urm(i)*(ury(i,2)-ury(i,1)));
 end
Rlax=0;
Rrax=0;
%% Calculation of support reaction at base of two column

% % Doing moment at left bottom of the column to be 0
% Rra = (-mar+mab+mal )/lb;
% 
% % Calculating load  Rla 
% Rla = fb - Rra;
% 
% % Calculating load Rlax
% Rlax = fr - fl;

if((s1==1 && s2==1)) %left fixed and right roller
    % Doing moment at left bottom of the column to be 0
  Rra = (mab+mal-marl)/lb;
  Rla = fb - Rra;% Calculating load  Rla 
   Rlax = fr - fl;% Calculating load Rla
   
   % Left Column
% Calculating moment at top of column
   Mlb =  -((Rlax*hl) + mbl) ; % A.C.W
Rlb = -Rla; %up positive
Rlbx = -fl-Rlax; % right positive
% right Column
Mrb = mbr; % ACW
Rrb  = -Rra; %up +ve
%Rrbx = fr;
elseif(s1==1 && s2==2)
    eqn1 = lb*Rra -mal +marl -mab -Rrax*(hl-hr) == 0;
    eqn2 = -lb*Rla +(hl-hr)*Rlax-malr +mbb + Mar == 0;
    eqn3 = Rra + Rla -fb == 0;
    eqn4 = Rla + Rrax +fl -fb ==0 ;
    sol = solve([eqn1, eqn2, eqn3, eqn4], [Rra, Rrax, Rla, Rlax]);

   % right Column
% Calculating moment at top of column
    Mlb =  -((Rlax*hl) + mbl) ; % A.C.W
    Rlb = -Rla; %up positive
    Rlbx = -fl-Rlax; % right positive
% right Column
    Mrb = mbr-Rrax*hl; % ACW
    Rrb  = -Rra; %up +ve
    Rrbx = fr;
end
%% Now dismantling frame into columns and beam

% % Left Column
% % Calculating moment at top of column
% Mlb =  Rlax*hl - mbl ;
% Rlb = Rla;
% Rlbx = -fl+Rlax;
% 
% % right Column
% Mrb = -mbr;
% Rrb  = -Rra;
% Rrbx = fr;

%% Calculating joint load on beam end
% Doing moment at left end is zero and calculating right end force

Rbb = abs(Rrb);
Rab = abs(Rlb);

%% Shear force vector for beam
     v = [Rab];
     cnt1=1;
     cnt2 = 1;
     cnt3=1;
     flag= 0 ;
     flag2 =0;
     flag3=0;
     for x = 0:0.02:lb
         if(cnt1<=npb && x==pbx(cnt1) && flag==0)
             v = [v v(end)-pbm(cnt1) ];
             cnt1 = cnt1 +1;
         end
         if(cnt2<=nub && (x>=ubx(cnt2,1)&& x<=ubx(cnt2,2)))

             %disp(xi);
             if(cnt3>1)
             v = [v v(end)-(ubm(cnt2)*((x-xi)))];
             flag3=1;
             else
                 flag=1;
             end
             cnt3 = cnt3 + 1;
         elseif(cnt2<=nub && x>ubx(cnt2,2))
                 cnt2 = cnt2 +1;
         else
           flag2=1;
           flag3=0;
                
         end
         if(cnt1<=npb && x==pbx(cnt1) && flag2==1)
                v = [v v(end)-pbm(cnt1) ];
                cnt1 = cnt1 +1;
                flag2=0;
         end
         if(x==lb)
             v = [v v(end)+Rbb];
         end
         %disp(x);
         %disp(v);
         if(flag3==0)
            v = [v v(end)];
         end
         xi = x;
     end
     v2=[];
     for x= 0:0.02:lb+0.02
         v2=[v2 0]; 
     end
     v1 = v(1:end-1);
     x= 0:0.02:lb+0.02;
     plot(x,v1,x,v2,'linewidth',2);
     xlabel('Length of the beam in m)');
     ylabel('Shear Force in KN');
     title('Shear force diagram');
     
 %% Shear force vector for left coulumn
     vl = [Rlax];
     cnt1=1;
     cnt2 = 1;
     cnt3=1;
     flag= 0 ;
     flag2 =0;
     flag3=0;
     for x = 0:0.1:hl
         if(cnt1<=npl && x==ply(cnt1) && flag==0)
             vl = [vl vl(end)-plm(cnt1) ];
             cnt1 = cnt1 +1;
         end
         if(cnt2<=nul && (x>=uly(cnt2,1)&& x<=uly(cnt2,2)))

             %disp(xi);
             if(cnt3>1)
             vl = [vl vl(end)-(ulm(cnt2)*((x-xi)))];
             flag3=1;
             else
                 flag=1;
             end
             cnt3 = cnt3 + 1;
         elseif(cnt2<=nul && x>uly(cnt2,2))
                 cnt2 = cnt2 +1;
         else
           flag2=1;
           flag3=0;
                
         end
         if(cnt1<=npl && x==ply(cnt1) && flag2==1)
                vl = [vl vl(end)-plm(cnt1) ];
                cnt1 = cnt1 +1;
                flag2=0;
         end
         if(x==hl)
             vl = [vl vl(end)-Rlbx];
         end
         disp(x);
         disp(v);
         if(flag3==0)
            vl = [vl vl(end)];
         end
         xi = x;
     end
     
     x= 0:0.1:hl+0.2;
     plot(x,vl);
     
 %% Shear force vector for right coulumn
     vr = [Rrax];
     cnt1=1;
     cnt2 = 1;
     cnt3=1;
     flag= 0 ;
     flag2 =0;
     flag3=0;
     for x = 0:0.1:hr
         if(cnt1<=npr && x==pry(cnt1) && flag==0)
             vr = [vr vr(end)-prm(cnt1) ];
             cnt1 = cnt1 +1;
         end
         if(cnt2<=nur && (x>=ury(cnt2,1)&& x<=ury(cnt2,2)))

             %disp(xi);
             if(cnt3>1)
             vr = [vr vr(end)-(urm(cnt2)*((x-xi)))];
             flag3=1;
             else
                 flag=1;
             end
             cnt3 = cnt3 + 1;
         elseif(cnt2<=nur && x>ury(cnt2,2))
                 cnt2 = cnt2 +1;
         else
           flag2=1;
           flag3=0;
                
         end
         if(cnt1<=npr && x==pry(cnt1) && flag2==1)
                vr = [vr vr(end)-prm(cnt1) ];
                cnt1 = cnt1 +1;
                flag2=0;
         end
         if(x==hr)
             vr = [vr vr(end)-Rrbx];
         end
         disp(x);
         disp(vr);
         if(flag3==0)
            vr = [vr vr(end)];
         end
         xi = x;
     end
     
     x= 0:0.1:hr+0.3;
     plot(x,vr);
     

%% shear force
sfd(10,1,0,10,5,5,5);
% j = 0;
% for x=0:0.02:lb
%     j = j + 1;
%     temp1 = 0;
%     for pi = 1:npb
%       if(x<pbx(pi))
%         v1(j) = 0;
%         break;
%       else
%           temp1 = temp1 + pbm(pi);
%         v1(j) = 0;
%         m1(j) = 0 ;
%       end
%     end
%     for ui = 1:nub
%       if(x<ubx(ui,1))
%         v1(j) = 0;
%         break;
%       elseif(x>=ubx(ui,1) && x<ubx(ui,2))
%         temp1 = temp1 + ubm(ui)*(x-ubx(ui,1));
%         v1(j) = 0;
%       elseif(x>=ubx(ui,2))
%           temp1 = temp1 + ubm(ui)*(ubx(ui,2) - ubx(ui,1));
%           v1(j) = 0;
%           m1(j) = 0 ;
%       end
%     end
%     if x<lb
%      v(j) = Rab - temp1;
%     else
%      v(j) = Rab + Rbb -temp1;
%     end
% end
% x=0:0.02:lb;
% subplot(3,1,1)
% plot(x,v,x,v1,'linewidth',2);
% title('shear force');
% xlabel('length of the beam in mm');
% ylabel('shear force in n');  
%      
%% Bending moment of beam
bmd(10,1,0,10,5,5,5,0);
% j = 0;
% for x=0:0.02:lb
%     j = j + 1;
%     temp = 0;
%     for pi = 1:npb
%       if(x<pbx(pi))
%         m1(j)= 0;
%         break;
%       else
%         temp = temp + (pbm(pi)*(x-pbx(pi)));
%         m1(j) = 0 ;
%       end
%     end
%     for ui = 1:nub
%       if(x<ubx(ui,1))
%         m1(j)= 0;   
%         break;
%       elseif(x>=ubx(ui,1) && x<ubx(ui,2))
%         temp = temp + ((ubm(ui)*((x-ubx(ui,1))^2))/2);
%         m1(j) = 0 ;
%       elseif(x>=ubx(ui,2))
%           temp = temp + (ubm(ui)*(ubx(ui,2) - ubx(ui,1))*(x-((ubx(ui,2) + ubx(ui,1))/2)));
%           m1(j) = 0 ;
%       end
%     end
%     m(j) = abs(Mlb)+(x*Rab) -temp;
% end
% x=0:0.02:lb;
% subplot(3,1,2)
% plot(x,m,x,m1,'linewidth',2);
% title('bending moment');
% xlabel('length of the beam in mm');
% ylabel('bending moment in n-mm');

%% Dividing udl for choosing element
% uc1 = [0,0,0];
% flag=0;
% i=1;
% while i<=nub-1
% disp(i);
%     if((ubx(i,2)<ubx(i+1,1))&&flag==0)
%         uc1=[uc1;ubx(i,1),ubx(i,2),ubm(i)];
%     elseif((ubx(i,2)>ubx(i+1,1))&&(ubx(i,2)<ubx(i+1,2)))
%         uc1 = [uc1; ubx(i,1),ubx(i+1,1),ubm(i);ubx(i+1,1),ubx(i,2),ubm(i)+ubm(i+1);ubx(i,2),ubx(i+1,2),ubm(i+1)];
%         %i=i+1;
%         disp('s');
%         flag=1;
%     elseif((ubx(i,2)>ubx(i+1,1))&&(ubx(i,2)>ubx(i+1,2)))
%         uc1 = [uc1; ubx(i,1),ubx(i+1,1),ubm(i);ubx(i+1,1),ubx(i+1,2),ubm(i)+ubm(i+1);ubx(i+1,2),ubx(i,2),ubm(i)];
%         flag=1;
%     elseif((ubx(i,2)<ubx(i+1,1))&&flag==1)
%         uc1=[uc1;ubx(i+1,1),ubx(i+1,2),ubm(i+1)];
%         flag=0;
%     end
%     i=i+1;
% end
% uc1 = [0,0,0];
% flag=0;
% i=1;
% while i<=nub-1
% disp(i);
%     if((ubx(i,2)<ubx(i+1,1)))
%         uc1=[uc1;ubx(i,1),ubx(i,2),ubm(i)];
%     elseif((ubx(i,2)>ubx(i+1,1))&&(ubx(i,2)<ubx(i+1,2)))
%         uc1 = [uc1; ubx(i,1),ubx(i+1,1),ubm(i);ubx(i+1,1),ubx(i,2),ubm(i)+ubm(i+1)];
%         %i=i+1;
%         disp('s');
% 
%     elseif((ubx(i,1)<ubx(i-1,2))&&(ubx(i,2)>ubx(i-1,2))&&i>1)
%         uc1 = [uc1; ubx(i-1,2),ubx(i,2),ubm(i)];
%         
%     elseif((ubx(i,2)>ubx(i+1,1))&&(ubx(i,2)>ubx(i+1,2)))
%         uc1 = [uc1; ubx(i,1),ubx(i+1,1),ubm(i);ubx(i+1,1),ubx(i+1,2),ubm(i)+ubm(i+1);ubx(i+1,2),ubx(i,2),ubm(i)];
% %         flag=1;
% %     elseif((ubx(i,2)<ubx(i+1,1))&&flag==1)
% %         uc1=[uc1;ubx(i+1,1),ubx(i+1,2),ubm(i+1)];
% %         flag=0;
%     end
%     i=i+1;
% end

uc1=[0,0,0];
i=1;
while(i<nub)
       if(i>1 && (ubx(i,1)>=ubx(i-1,2) && ubx(i,2)>ubx(i-1,2))&&(ubx(i,2)<=ubx(i+1,1))&&(ubx(i,2)<ubx(i+1,2)))
      uc1=[uc1;ubx(i,1),ubx(i,2),ubm(i)];
       elseif(i<=1 &&(ubx(i,2)<=ubx(i+1,1))&&(ubx(i,2)<ubx(i+1,2)))
           uc1=[uc1;ubx(i,1),ubx(i,2),ubm(i)];
       
    elseif((ubx(i,2)>ubx(i+1,1))&&(ubx(i,2)<ubx(i+1,2)))
      uc1 = [uc1; ubx(i,1),ubx(i+1,1),ubm(i);ubx(i+1,1),ubx(i,2),ubm(i)+ubm(i+1)];  
    elseif(i>1&&(ubx(i,1)<ubx(i-1,2))&&(ubx(i,2)>ubx(i-1,2)))
        uc1 = [uc1; ubx(i-1,2),ubx(i,2),ubm(i)];
    elseif((ubx(i,1)<=ubx(i+1,1))&&(ubx(i,2)>=ubx(i+1,2)))
        uc1 = [uc1; ubx(i,1),ubx(i+1,1),ubm(i);ubx(i+1,1),ubx(i+1,2),ubm(i)+ubm(i+1);ubx(i+1,2),ubx(i,2),ubm(i)];
        i=i+1;
    end
    i=i+1;
end
%for adding udl if it is one
if(i==nub)
      uc1=[uc1;ubx(i,1),ubx(i,2),ubm(i)];
end

if( i>1 &&(ubx(i,1)<ubx(i-1,2))&&(ubx(i,2)>ubx(i-1,2)))
        uc1 = [uc1; ubx(i-1,2),ubx(i,2),ubm(i)];
elseif( i>1 &&(ubx(i-1,2)<=ubx(i,1))&&(ubx(i-1,2)<ubx(i,2)))
    uc1=[uc1;ubx(i,1),ubx(i,2),ubm(i)];
end
for i = 1:(size(uc1,1)-1)
    if(uc1(i,2)<uc1(i+1,1))
        uc1=[uc1;uc1(i,2),uc1(i+1,1),0];
    end
end
if(uc1(end,2)<lb &&(ubx(end,2)<uc1(end,2)))
    uc1=[uc1;uc1(end,2),lb,0];
end

uc1 = uc1(2:end,:);
uc1 = sortrows(uc1,1);

% for dividing element according to the point load 
j=1;
uc=[0,0,0,0];
for i=1:size(uc1)
    if(j<=npb&&(uc1(i,1)<pbx(j))&&(pbx(j)<uc1(i,2)))
        disp(pbx(j));
        disp(uc1(i,1));
        uc=[uc;uc1(i,1),pbx(j),uc1(i,3),pbm(j);pbx(j),uc1(i,2),uc1(i,3),pbm(j)];
        j=j+1;
    else
        uc=[uc;uc1(i,1),uc1(i,2),uc1(i,3),0];
    end
end



% j=1;
% k=zeros(4,4,size(uc1,1));
% f = zeros(4,1,size(uc1,1));
% for i  = 1:size(uc1,1)
%     if(uc1(i,1)<=pbx(j) && uc1(i,2)<pbx(j))
%         if(uc1(i,3) == 0 )
%            k(:,:,i)=((E*I)/(uc1(i,2)-uc1(i,1))^3)*[12,6*(uc1(i,2)-uc1(i,1)),-12,6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),4*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2);-12,-6*(uc1(i,2)-uc1(i,1)),12,-6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),4*(uc1(i,2)-uc1(i,1))]; 
%            f(:,1,i)=[0;0;0;0];
%         else
%            k(:,:,i)=((E*I)/(uc1(i,2)-uc1(i,1))^3)*[12,6*(uc1(i,2)-uc1(i,1)),-12,6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),4*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2);-12,-6*(uc1(i,2)-uc1(i,1)),12,-6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),4*(uc1(i,2)-uc1(i,1))]; 
%            f(:,1,i)=[uc1(i,3)*((uc1(i,2)-uc1(i,1))*0.5);(-uc1(i,3)/12)*((uc1(i,2)-uc1(i,1))^2);uc1(i,3)*((uc1(i,2)-uc1(i,1))*0.5);(uc1(i,3)/12)*((uc1(i,2)-uc1(i,1))^2)];
%         end
%     if(uc1(i,1)<pbx(j) && uc1(i,2)>=pbx(j))
%         a = pbx(j) - uc1(i,1);
%         b = uc1(i,2) - pbx(j);
%         l = uc1(i,2)-uc1(i,1);
%         if(uc1(i,3)==0)
%            k(:,:,i)=((E*I)/(uc1(i,2)-uc1(i,1))^3)*[12,6*(uc1(i,2)-uc1(i,1)),-12,6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),4*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2);-12,-6*(uc1(i,2)-uc1(i,1)),12,-6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),4*(uc1(i,2)-uc1(i,1))]; 
%            f(:,1,i)=[((pbm(j)*b*b)/(l^3))*((3*a)+b);(pbm(j)*a*b*b)/(l^2);((pbm(j)*a*a)/(l^3))*((3*b)+a);(pbm(j)*a*a*b)/(l^2)];
%         else
%            k(:,:,i)=((E*I)/(uc1(i,2)-uc1(i,1))^3)*[12,6*(uc1(i,2)-uc1(i,1)),-12,6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),4*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2);-12,-6*(uc1(i,2)-uc1(i,1)),12,-6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),4*(uc1(i,2)-uc1(i,1))]; 
%            f(:,1,i)=[(uc1(i,3)*((uc1(i,2)-uc1(i,1))*0.5))+(((pbm(j)*b*b)/(l^3))*((3*a)+b));((-uc1(i,3)/12)*((uc1(i,2)-uc1(i,1))^2))+((pbm(j)*a*b*b)/(l^2));(uc1(i,3)*((uc1(i,2)-uc1(i,1))*0.5))+(((pbm(j)*a*a)/(l^3))*((3*b)+a));((uc1(i,3)/12)*((uc1(i,2)-uc1(i,1))^2))+((pbm(j)*a*a*b)/(l^2))];
%         end
%         if(j<npb)
%             j=j+1;
%         end
%     else
%           if(uc1(i,3) == 0 )
%            k(:,:,i)=((E*I)/(uc1(i,2)-uc1(i,1))^3)*[12,6*(uc1(i,2)-uc1(i,1)),-12,6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),4*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2);-12,-6*(uc1(i,2)-uc1(i,1)),12,-6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),4*(uc1(i,2)-uc1(i,1))]; 
%            f(:,1,i)=[0;0;0;0];
%         else
%            k(:,:,i)=((E*I)/(uc1(i,2)-uc1(i,1))^3)*[12,6*(uc1(i,2)-uc1(i,1)),-12,6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),4*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2);-12,-6*(uc1(i,2)-uc1(i,1)),12,-6*(uc1(i,2)-uc1(i,1));6*(uc1(i,2)-uc1(i,1)),2*((uc1(i,2)-uc1(i,1))^2),-6*(uc1(i,2)-uc1(i,1)),4*(uc1(i,2)-uc1(i,1))]; 
%            f(:,1,i)=[uc1(i,3)*((uc1(i,2)-uc1(i,1))*0.5);(-uc1(i,3)/12)*((uc1(i,2)-uc1(i,1))^2);uc1(i,3)*((uc1(i,2)-uc1(i,1))*0.5);(uc1(i,3)/12)*((uc1(i,2)-uc1(i,1))^2)];
%         end
%     end
% end
% 
% % F=KQ
% syms q1 q2 q3 q4 Q
% for i=1:size(uc1,1)
%     l = uc1(i,2) - uc1(i,1);
%     v1 = ((6*E*I)/(l^3))*((2*q1) + (l*q2)-(2*q3)+(l*q4));
%     m1 = ((E*I)/(l^2))*((-6*q1)+(-4*l*q2) +(6*q3)-(2*l*q4));
%     v2 = ((6*E*I)/(l^3))*((2*q1) + (l*q2)-(2*q3)+(l*q4));
%     m2 = ((E*I)/(l^2))*((6*q1)+(2*l*q2) +(-6*q3)+(4*l*q4));    
%     eq1 =   (k(1,:,i)*[q1;q2;q3;q4]) + f(1,1,i) - v1 == 0; 
%     eq2 =   (k(2,:,i)*[q1;q2;q3;q4]) + f(2,1,i) - m1 == 0;
%     eq3 =   (k(3,:,i)*[q1;q2;q3;q4]) + f(3,1,i) - v2 == 0;
%     eq4 =   (k(4,:,i)*[q1;q2;q3;q4]) + f(4,1,i) - m2 == 0;
%     sol = solve([eq1, eq2, eq3, eq4], [q1, q2, q3, q4]);
%     Q = [Q; sol.q1,sol.q2,sol.q3,sol.q4];
% end

