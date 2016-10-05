%% shear force
function sfd(lb,npb,nub,pbm,pbx,Rab,Rbb)
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
assignin('base','v',v);
x=0:0.02:lb;
subplot(3,1,1)
plot(x,v,x,v1,'linewidth',2);
title('shear force');
xlabel('length of the beam in mm');
ylabel('shear force in n'); 
end