%%bending moment;
function bmd(lb,npb,nub,pbm,pbx,Rab,Rbb,Mlb)
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
assignin('base','m',m);
x=0:0.02:lb;
subplot(3,1,2)
plot(x,m,x,m1,'linewidth',2);
title('bending moment');
xlabel('length of the beam in mm');
ylabel('bending moment in n-mm');
