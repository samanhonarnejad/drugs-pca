function crosshairs(pos1, pos2)
mn = 0.5 * (pos1 + pos2);
line([pos1(1), pos2(1)], [mn(2), mn(2)], 'linewidth', 1, 'color', 'k');
line([mn(1), mn(1)], [pos1(2), pos2(2)], 'linewidth', 1, 'color', 'k');
