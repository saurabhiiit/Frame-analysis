
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
 Rbb = abs(mab)/lb;
 Rab = abs(fb) - Rbb;

    %% Shear force vector
     v = [Rab];
     cnt1=1;
     cnt2 = 1;
     cnt3=1;
     flag= 0 ;
     flag2 =0;
     flag3=0;
     for x = 0:0.1:lb
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
         disp(x);
         disp(v);
         if(flag3==0)
            v = [v v(end)];
         end
         xi = x;
     end
     
     x= 0:0.1:lb+0.3;
     plot(x,v);