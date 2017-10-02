function ys=locsmooth(y,span)
% function that performs local smoothing (moving window average)

  if (mod (span,2) == 0)
    error ('smooth: span must be odd.')
  end
  ys = y;
  l = length(y);
  for i=1:l
    a=i-(span-1)/2;
    b=i+(span-1)/2;
    if a>0 && b<=l
        ys(i)=mean(y(a:b));
    elseif a<=0
		ys(i)=mean(y(1:b));
    else
	    ys(i)=mean(y(a:l));
    end
  end
