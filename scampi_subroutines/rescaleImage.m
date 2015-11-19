function [rescaledImage] = rescaleImage(original_image)

rescaledImage = (original_image - min(original_image(:) ) ) ./ (max(original_image(:) ) - min(original_image(:) ) );

end