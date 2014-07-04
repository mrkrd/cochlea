function resampled = resample(X, P, Q)

  warning('Using fake resample().')

  b = fir1(Q, 1/Q);
  a = [1];

  if (Q > 1) & (P == 1)
    filtered = filtfilt(b, a, X);
    resampled = filtered(1:Q:end);
  elseif (Q == 1) & (P == 1)
    resampled = X;
  else
    error('resample: unnown case')
  end


  %else
  %  resampled = repmat(X, P, 1);
  %  resampled = resampled(:)';
  %end
