using Images
function correlate_image(img, template)
    z = imfilter(Float64, img, centered(template), Algorithm.FIR())
    q = sqrt.(imfilter(Float64, img.^2, centered(ones(size(template))),
                       Algorithm.FIR()))
    m = sqrt(sum(template.^2))
    corr = z ./ (m * max.(eps(), q))
    return corr
end
