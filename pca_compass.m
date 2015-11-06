function pca_compass(coeff, names, arr_zero, arr_len, mag)
idx_stain = find(sqrt(sum(coeff(:, [1, 2]) .^ 2, 2)) > mag);
for sig = idx_stain'
    vec = zeros(size(coeff, 1), 1);
    vec(sig) = arr_len;
    vec = vec' * coeff;
    h = annotation('arrow');
    set(h, 'parent', gca(), 'position', [arr_zero, vec(1), vec(2)]);
    text(1.2 * vec(1) + arr_zero(1), 1.2 * vec(2) + arr_zero(2), ...
        names{sig});
end
