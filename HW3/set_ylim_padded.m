function set_ylim_padded(data, pad_frac)
    if nargin < 2
        pad_frac = 0.1;  % 10% of the data's range, a reasonable default
    end
    lo = min(data);
    hi = max(data);
    range = hi - lo;
    if range == 0
        % constant signal (e.g. all zeros) — pad by an absolute amount instead
        pad = max(abs(hi), 1) * 0.1;
    else
        pad = range * pad_frac;
    end
    ylim([lo - pad, hi + pad]);
end