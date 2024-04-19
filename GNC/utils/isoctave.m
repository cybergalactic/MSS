function is_oct = isoctave()
    is_oct = exist('OCTAVE_VERSION', 'builtin') > 0;
end