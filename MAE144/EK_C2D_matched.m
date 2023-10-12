function [bz, az, K] = EK_C2D_matched(bs,as,h,omega, proper)
%[bz, az, K] = EK_C2D_matched(bs,as,h,omega, proper)
%Convert D(s)=bs(s)/as(s) to D(z)=K*bz(z)/az(z) using matched z-transform approach.
%
%Returns the coefficient vectors for bz and az, as well as the matched gain
%(K) for the transformation.
%
%bs and as should be inputs of coefficient vectors for D(s). Symbolic
%coefficients will also be accepted.
%
%h is the timestep we are looking at.
%
%omega is an optional input. If inputted it will do the transformation and
%match the gains at the inputted omega. Otherwise omega is assumed to be 0.
%
%proper is an optional input. It takes a true or false value. If true is
%entered, the transformation will be done and return a strictly causal
%D(z). If false is enetered, a semi-causal D(z) will also be an acceptable
%output of the function. 
%
%A simple test to show that this function obtains the same answer as c2d
%with the 'matched' option:
%
%>> [bz, az, K] = EK_C2D_matched([1 1], [1 10 0], 1, 0, true)
%>> c2d(tf([1 1], [1 10 0]), 1, 'matched')
%Running these will give the same D(z)!

    
    %Check if we are given an omega of interest, otherwise just set it to 0
    if ~exist('omega','var')
      omega = 0.0000;
    end

    %Check if we are told whether we are ok with semi-causal transfer
    %functions. If not, we default to being ok with semi-causality.
    if ~exist('proper','var')
      proper = false;
    end

    %Find the poles of as
    as_poles = RR_roots(RR_poly(as));
    %Find the zeros of bz by using our conversion on the zeros of bs
    bz_zeros = exp(RR_roots(RR_poly(bs))*h);
    %Find the poles of az by using our conversion on the poles of as
    az_poles = exp(RR_roots(RR_poly(as))*h);
    
    %If we end up with any super small values in the roots of az and bz,
    %set them to 0 since this is likely due to rounding. 
    bz_zeros(isAlways(abs(bz_zeros)<=0.0001)) = 0;
    az_poles(isAlways(abs(az_poles)<=0.0001)) = 0;

    %Check to see if either of our denominators will evaluate to 0 when
    %trying to find the gain. If it will, add 0.1 to omega so we are still
    %in the vicinity of the omega of interest, but will have a finite gain.
    while ismember(exp(omega*h*j),az_poles) || ismember(omega*j,as_poles)
        omega = omega + 0.1;
    end
    
    %Find difference in the order of the denominator and numerator of D(s).
    n = length(as) - length(bs);

    %Map infiite zeros of D(s) to z = -1 in D(z).
    for i = 0:n-1
        bz_zeros(end+1) = -1;
    end

    %If we want strict causality, and we added zeros at z = -1, remove one
    %of them.
    if proper == true && n > 0
        bz_zeros = bz_zeros(1:end-1);
    end
    
    %Extract the polynomial coefficient vectors fo bz and az from our
    %transofrmed poles and zeros. 
    bz = RR_poly(bz_zeros,1).poly;
    az = RR_poly(az_poles,1).poly;
    
   
    %This line helps to simplify the numbers when we have symbolics.
    sympref('FloatingPointOutput', true);

    %Evaluate the gains of D(s) and D(z) at our omega value.
    Ds = abs(RR_evaluate(RR_tf(bs,as),j*omega*h));
    Dz = abs(RR_evaluate(RR_tf(bz_zeros,az_poles,1), exp(j*omega*h)));

    %This coefficient matches the gain of D(z) to D(s) at omega. 
    K = Ds/Dz;
    
end

