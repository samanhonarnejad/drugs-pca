function ellipse2d(pos1, pos2, col, alpha)
mn = 0.5 * (pos1 + pos2);
rx = pos1(1) - mn(1);
ry = pos1(2) - mn(2);
rad = 0 : (2 * pi / 32) : 2 * pi;
x = mn(1) + rx * sin(rad);
y = mn(2) + ry * cos(rad);
patch(x, y, col, 'facealpha', alpha, 'linestyle', 'none');
