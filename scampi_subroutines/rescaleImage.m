function [rescaledImage] = rescaleImage(original_image)

rescaledImage = (original_image - min(original_image(:) ) ) ./ (max(original_image(:) ) - min(original_image(:) ) );

rescaledImage2 = original_image - min(original_image(:) ); 
rescaledImage2 = rescaledImage2 ./ max(rescaledImage2(:) );

if sum(rescaledImage(:) ~= rescaledImage2(:) ) > 0; error; end

end