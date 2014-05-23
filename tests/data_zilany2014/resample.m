function resampled = resample_fake(X, P, Q)

  warning('Using fake resample().')

  b = fir1(Q, 1/Q);
  a = [1];

  filtered = filtfilt(b, a, X);
  resampled = filtered(1:Q:end);
