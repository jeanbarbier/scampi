function [dumped] = dumping(old, new, m)

% compute the dumping of b and c by a factor m
dumped = m .* old + (1 - m) .* new;

end