% Parameter Estimation using Recursive least squares
% Parallel predictor

clear

% actual values of the parameters
disp('**** PARALLEL PREDICTOR ****');
disp('         b0*z+b1');
disp('   G(z)=-----------');
disp('        z^2+a1*z+a2');
disp('input the actual parameter values and hit return');

% a1=input('a1=');
% a2=input('a2=');
% b0=input('b0=');
% b1=input('b1=');

a1=0.7;
a2=0.2;
b0=2;
b1=0.9;

disp('');
disp('**** select the input function ****');
disp('if u=random, then type -1');
disp('if u=sin(w1*k)+sin(w2*k)+sin(w3*k), then type 0');
disp('if u=1, then type 1');

% ct=input('  :  ');
ct=1;
if ct==0,
  disp('**** input values of w1, w2, w3 ****');
  w1=2;
  w2=5;
  w3=7;
end
ct1=0;
% while ct1==1,
%   disp('**** input forgetting factor ****');
%   disp('**** lamda1=1, lamda2=0: constant');  
%   disp('**** lamda1=1, lamda2=1: least squares');
%   disp('**** lamda1<1, lamda2=1: weighted least squares');
%   lam1=input('lamda1=');
%   lam2=input('lamda2=');
  lam1=1;
  lam2=0;

  % initial gain matrix
  disp('**** input the diagonal elements of the gain matrix ****');
  disp('**** F=diag(f11, f22, f33, f44) ****');
  fratio=1;
  fratio2=1;
%   f11=input('f11=');
%   f22=input('f22=');
%   f33=input('f33=');
%   f44=input('f44=');
  f11=fratio;
  f22=fratio;
  f33=fratio2;
  f44=fratio2;
%   
  F=diag([f11, f22, f33, f44]);
  disp('**** input gain ****');
  
  % estimates of the parameters
  thhat=[0.8 0.8 0.8 0.8]';
%     thhat=[0.8 0.8 0.8 0.8]';
  
  % initial data values
  % y1=y(k-1), y2=y(k-2), u1=u(k-1), u2=u(k-2)
  y1=0;y2=0;u1=0;u2=0;y1h=0;y2h=0;
  phi=[-y1h -y2h u1 u2]';
  
%   disp('');
%   disp('**** turn on or off the white noise ****');
%   disp('if you want to turn on, then type 1');
%   disp('if you want to turn off, then type 0');
%   ct2=input('  :  ');
  ct2=0;
  kn=600;       % final time
  if ct2==1,           
    wn1=0; wn0=0; %initial values for white noise effect w(0) and w(-1)
    wn=0.1*randn(kn);
  end
  
  for k=1:kn,
    % make new measurements
    if ct==0,
        inp='sinusoidal';
      u=sin(w1*k)+sin(w2*k)+sin(w3*k);
    elseif ct==1,
        inp='1';
      u=1;
    elseif ct==-1,
      u=randn;
    end

    if ct2==1,
       wn2=wn(k);
       y=-a1*y1-a2*y2+b0*u1+b1*u2+wn2+a1*wn1+a2*wn0;
       wn0=wn1;
       wn1=wn2;
   else
       y=-a1*y1-a2*y2+b0*u1+b1*u2;
   end
    
    %a priori predicted output and prediction error
    yh0=thhat'*phi;
    e0=y-yh0;
    
    % update parameter values
    gamma=F*phi/(1+phi'*F*phi);
    thhat=thhat+gamma*e0;
    
    % a posteriori predicted output and prediction error
    yh=thhat'*phi;
    e=y-yh;
    
    % update the gain matrix
    if lam2>0,
       F=(F-F*phi*phi'*F/(lam1/lam2+phi'*F*phi))/lam1;
    end

    
    % prepare for next recursions  
    y2h=y1h;
    y1h=yh;
    y2=y1;
    y1=y;
    u2=u1;
    u1=u;
    phi=[-y1h -y2h u1 u2]';
    
    % store values for plotting
    % thetahat: parameter estimate vector   
    thetahat(k, 1)=thhat(1); thetahat(k, 2)=thhat(2);
    thetahat(k, 3)=thhat(3); thetahat(k, 4)=thhat(4);

    yr(k)=y;        % real output y
    yhat0(k)=yh0;   % a priori predicted output
    err0(k)=e0;     % a priori prediction error
    yhat(k)=yh;     % a posteriori predicted output
    err(k)=e;       % a posteriori prediction error
    f(k)=det(F);    % determinant of gain matrix
    t(k)=k;         % time
 
  end
  
  % plotting error, gain matrix and estimates
  figure;
  subplot(211)
  plot(t,err); title('prediction error, parallel predictor');
  subplot(212)
  plot(t, f); title('determinant of gain matrix');
  if ct==0,
    xlabel('input: sum of sinusoids');
  elseif ct==1, 
    xlabel('input: constant');
  elseif ct==-1,
    xlabel('input:random');
    
  end % ending if
  
  disp('pause ...'); pause
  figure;
  subplot(221); plot(t, thetahat(:, 1)); ylabel('a_1');
  subplot(222); plot(t, thetahat(:, 2)); ylabel('a_2');
  subplot(223); plot(t, thetahat(:, 3)); ylabel('b_0');
  subplot(224); plot(t, thetahat(:, 4)); ylabel('b_1');
  a1_hat=thetahat(length(t),1)
  a2_hat=thetahat(length(t),2)
  b0_hat=thetahat(length(t),3)
  b1_hat=thetahat(length(t),4) 
  theta_hat=[a1_hat a2_hat b0_hat b1_hat]
  sgtitle(["input="+inp+","+"lamda1="+lam1+","+"lamda2="+lam2+","+"f ratio="+fratio+"/"+fratio2,"theta_hat=["+theta_hat(1,1)+","+theta_hat(1,2)+","+theta_hat(1,3)+","+theta_hat(1,4)+"]"]);
%   ct1=input('do you wish to try different lamda values, type 1 else 0: ');
%   ct1=0;
% end






