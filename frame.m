%%% Design of frame
clc

%% input for beam 
lb = input('Enter the span of the beam: ');

npb = input('Enter the number of the point load acting: ');
mab = 0;
fb = 0;

for i =1:npb
    pbx(i) = input('Enter the coodinate where point load is acting ');
    pbm(i) = input('Enter the point load magnitute on beam: ');
    mab = mab + (pbx(i)*pbm(i));
    fb = fb + pbm(i);
end

nub = input('Enter the number of UDL acting on the beam: ');
for i = 1:nub
    ubx(i,1) = input('Enter the left end of the UDL coodinate where point load is acting: ');
    ubx(i,2) = input('Enter the right end of the UDL coodinate where point load is acting: ');
    ubm(i) = input('Enter the UDL acting on the beam: ');
    mab = mab + (ubm(i)*(ubx(i,2)-ubx(i,1))*0.5*(ubx(i,2)+ubx(i,1)));
    fb = fb + (ubm(i)*(ubx(i,2)-ubx(i,1)));
end

%% Input for left column
disp('Enter the detail of left end of the coulumn:-');
hl = input('Enter the height of the left end of coulumn: ');    

mal = 0;
fl = 0 ;
mbl=0;
npl = input('Enter the number of point load in left coulumn: ');
for i = 1:npl
    plm(i) = input('Enter the magnitude of point load acting: ');
    ply(i) = input('Enter the coordinate where point load is acting: ');
    mal = mal + (plm(i)*ply(i));
    mbl = mbl + (plm(i)*(hl-ply(i)));
    fl = fl + plm(i);
end

nul = input('Enter the number of UDL load in left coulumn: ');
for i = 1:nul
    ulm(i) = input('Enter the magnitude of UDL load acting: ');
    uly(i,1) = input('Enter the bottom coordinate where UDL load is started acting: ');
    uly(i,2) = input('Enter the upper coordinate where UDL load is ended acting: ');    
    mal = mal + (ulm(i)*(uly(i,2)-uly(i,1))*0.5*(uly(i,2)+uly(i,1)));
    mbl = mbl + (ulm(i)*(uly(i,2)-uly(i,1))*0.5*(hl-(uly(i,2)+uly(i,1))));
    fl = fl + (ulm(i)*(uly(i,2)-uly(i,1)));
end

%% Input for right column
disp('Enter the detail of right end of the coulumn:-');
hr = input('Enter the height of the right end of coulumn: ');    

 mar = 0;
 fr = 0 ;
 mbr=0;
 npr = input('Enter the number of point load in right coulumn: ');
 for i = 1:npr
     prm(i) = input('Enter the magnitude of point load acting: ');
     pry(i) = input('Enter the coordinate where point load is acting: ');
     mar = mar + (prm(i)*pry(i));
     mbl = mbr + (prm(i)*(hr-pry(i)));
     fr = fr + prm(i);
 end
 
 nur = input('Enter the number of UDL load in left coulumn: ');
 for i = 1:nur
     urm(i) = input('Enter the magnitude of UDL load acting: ');
     ury(i,1) = input('Enter the bottom coordinate where UDL load is started acting: ');
     ury(i,2) = input('Enter the upper coordinate where UDL load is ended acting: ');    
     mar = mar + (urm(i)*(ury(i,2)-ury(i,1))*0.5*(ury(i,2)+ury(i,1)));
     mbr = mbr + (urm(i)*(ury(i,2)-ury(i,1))*0.5*(hr-(ury(i,2)+ury(i,1))));
     fr = fr + (urm(i)*(ury(i,2)-ury(i,1)));
 end

%% Calculation of support reaction at base of two column

% Doing moment at left bottom of the column to be 0
Rra = (-mar+mab+mal )/lb;

% Calculating load  Rla 
Rla = fb - Rra;

% Calculating load Rlax
Rlax = fr - fl;

%% Now dismantling frame into columns and beam

% Left Column
% Calculating moment at top of column
Mlb =  Rlax*hl - mbl ;
Rlb = Rla;
Rlbx = -fl+Rlax;

% right Column
Mrb = -mbr;
Rrb  = -Rra;
Rrbx = fr;

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
j = 0;
for x=0:0.02:lb
    j = j + 1;
    temp1 = 0;
    for pi = 1:npb
      if(x<pbx(pi))
        v1(j) = 0;
        break;
      else
          temp1 = temp1 + pbm(pi);
        v1(j) = 0;
        m1(j) = 0 ;
      end
    end
    for ui = 1:nub
      if(x<ubx(ui,1))
        v1(j) = 0;
        break;
      elseif(x>=ubx(ui,1) && x<ubx(ui,2))
        temp1 = temp1 + ubm(ui)*(x-ubx(ui,1));
        v1(j) = 0;
      elseif(x>=ubx(ui,2))
          temp1 = temp1 + ubm(ui)*(ubx(ui,2) - ubx(ui,1));
          v1(j) = 0;
          m1(j) = 0 ;
      end
    end
    if x<lb
     v(j) = Rab - temp1;
    else
     v(j) = Rab + Rbb -temp1;
    end
end
x=0:0.02:lb;
subplot(3,1,1)
plot(x,v,x,v1,'linewidth',2);
title('shear force');
xlabel('length of the beam in mm');
ylabel('shear force in n');  
     
%% Bending moment of beam
j = 0;
for x=0:0.02:lb
    j = j + 1;
    temp = 0;
    for pi = 1:npb
      if(x<pbx(pi))
        m1(j)= 0;
        break;
      else
        temp = temp + (pbm(pi)*(x-pbx(pi)));
        m1(j) = 0 ;
      end
    end
    for ui = 1:nub
      if(x<ubx(ui,1))
        m1(j)= 0;   
        break;
      elseif(x>=ubx(ui,1) && x<ubx(ui,2))
        temp = temp + ((ubm(ui)*((x-ubx(ui,1))^2))/2);
        m1(j) = 0 ;
      elseif(x>=ubx(ui,2))
          temp = temp + (ubm(ui)*(ubx(ui,2) - ubx(ui,1))*(x-((ubx(ui,2) + ubx(ui,1))/2)));
          m1(j) = 0 ;
      end
    end
    m(j) = abs(Mlb)+(x*Rab) -temp;
end
x=0:0.02:lb;
subplot(3,1,2)
plot(x,m,x,m1,'linewidth',2);
title('bending moment');
xlabel('length of the beam in mm');
ylabel('bending moment in n-mm');
