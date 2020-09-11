function g = deconv_conv(f,t,delta)

  if nargin==2,delta=false;end

  T=0.1; sig=0.06;

  dt = mean(diff(t));
  %dt = t(2)-t(1); 
  nt = length(t);

  % convolve analytical solution (for delta source) with STF
  if delta
    filter=(1/T^2)*t.*exp(-t./T);
    f=dt*conv(f,filter);
    f=f(1:nt);
  end  

  ts=4*sig; tau=t-ts;
  factor=1-(2*T/sig^2)*tau-((T/sig)^2)*(1-(tau./sig).^2);
  filter=(1/sqrt(2*pi)/sig)*factor.*exp(-0.5*(tau./sig).^2);
  g=dt*conv(f,filter);
  g=g(1:nt);

