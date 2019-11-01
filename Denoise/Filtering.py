import numpy as np
import math



class Filtering:
    image = None
    filter_name = None
    filter_size = None
    alpha_d = None
    order = None
    np.seterr(divide='ignore', invalid='ignore')

    def __init__(self, image, filter_name, filter_size, alpha_d=2, order = 1.5 ):
        """initializes the variables of spatial filtering on an input image
        takes as input:
        image: the noisy input image
        filter_name: the name of the mask to use
        filter_size: integer value of the size of the mask
        alpha_d: parameter of the alpha trimmed mean filter
        order: parameter of the order for contra harmonic"""

        self.image = image
        if filter_name == 'median':
            self.filter = self.get_median_filter
        elif filter_name == 'min':
            self.filter = self.get_min_filter
        elif filter_name == 'max':
            self.filter = self.get_max_filter
        elif filter_name == 'alpha_trimmed':
            self.filter = self.get_alpha_filter
        elif filter_name == 'arithmetic_mean':
            self.filter = self.get_arithmetic_mean_filter
        elif filter_name == 'geometric_mean':
            self.filter = self.get_geo_mean_filter
        elif filter_name == 'contra_harmonic':
            self.filter = self.get_contra_harmonic_filter

        self.filter_size = filter_size
        self.alpha_d = alpha_d
        self.order = order
        self.filter_name = filter_name


    def get_median_filter(self, kernel):
        """Computes the median filter
        takes as input:
        kernel: a list/array of intensity values
        returns the median value in the current kernel
        """

        return np.median(kernel)

        #return 0

    def get_min_filter(self, kernel):
        """Computes the minimum filter
        takes as input:
        ikernel: a list/array of intensity values
        returns the minimum value in the current kernel"""

        return np.min(kernel)

        #return 0

    def get_max_filter(self, kernel):
        """Computes the maximum filter
        takes as input:
        kernel: a list/array of intensity values
        returns the maximum value in the current kernel"""

        return np.max(kernel)

        #return 0

    def get_alpha_filter(self, kernel, alpha_d):
        """Computes the median filter
        takes as input:
        kernel: a list/array of intensity values
        alpha_d: clip off parameter for the alpha trimmed filter
        returns the alpha trimmed mean value in the current kernel"""

        for i in range(0, math.floor(alpha_d/2)):
            kernel.remove(min(kernel))
            kernel.remove(max(kernel))
        return np.mean(kernel)

        #return 0

    def get_arithmetic_mean_filter(self, kernel):
        """Computes the arithmetic mean filter
        takes as input:
        kernel: a list/array of intensity values
        returns the arithmetic mean value in the current kernel"""

        return np.mean(kernel)

        #return 0

    def get_geo_mean_filter(self, kernel):
        """Computes the geometric mean filter
                        takes as input:
                        kernel: a list/array of intensity values
                        returns the geometric mean value in the current kernel"""

        geo_mean = math.pow(np.prod(kernel), 1/len(kernel))
        return geo_mean

        #return 0

    def get_contra_harmonic_filter(self, kernel, order):
        """Computes the harmonic filter
                        takes as input:
        kernel: a list/array of intensity values
        order: order paramter for the
        returns the harmonic mean value in the current kernel"""
        har_mean = sum(np.power(kernel, order+1)) / sum(np.power(kernel, order))
        return har_mean

        #return 0

    def filtering(self):
        """performs filtering on an image containing gaussian or salt & pepper noise
        returns the denoised image
        ----------------------------------------------------------
        Note: Here when we perform filtering we are not doing convolution.
        For every pixel in the image, we select a neighborhood of values defined by the kernal and apply a mathematical
        operation for all the elements with in the kernel. For example, mean, median and etc.

        Steps:
        1. add the necesssary zero padding to the noisy image, that way we have sufficient values to perform the operati
        ons on the pixels at the image corners. The number of rows and columns of zero padding is defined by the kernel size
        2. Iterate through the image and every pixel (i,j) gather the neighbors defined by the kernel into a list (or any data structure)
        3. Pass these values to one of the filters that will compute the necessary mathematical operations (mean, median, etc.)
        4. Save the results at (i,j) in the ouput image.
        5. return the output image
        """

        r,c = self.image.shape
        size_1 = self.filter_size
        size_1 = np.int(0.5 * size_1 - 0.5)
        pad_image = np.zeros((r + (2 * size_1), c + (2 * size_1)))
        ker = []

        for i in range(c):
            for j in range(r):
                pad_image[i+size_1, j+size_1] = self.image[i,j]

        output_img = np.zeros([r, c])
        for i in range(r):
            for j in range(c):
                new_i = i+1
                new_j = j+1
                ker.clear()
                for x in range(-size_1,size_1 + 1):
                    for y in range(-size_1,size_1 + 1):
                        ker.append(pad_image[new_i+x, new_j+y])
                if self.filter_name == 'alpha_trimmed':
                    mask = self.filter(ker, self.alpha_d)
                elif self.filter_name == 'contra_harmonic':
                    mask = self.filter(ker, self.order)
                else:
                    mask = self.filter(ker)
                output_img[i][j] = mask

        return output_img

# def filtering(self):
#         """performs filtering on an image containing gaussian or salt & pepper noise
#         returns the denoised image
#         ----------------------------------------------------------
#         Note: Here when we perform filtering we are not doing convolution.
#         For every pixel in the image, we select a neighborhood of values defined by the kernal and apply a mathematical
#         operation for all the elements with in the kernel. For example, mean, median and etc.
#
#         Steps:
#         1. add the necesssary zero padding to the noisy image, that way we have sufficient values to perform the operati
#         ons on the pixels at the image corners. The number of rows and columns of zero padding is defined by the kernel size
#         2. Iterate through the image and every pixel (i,j) gather the neighbors defined by the kernel into a list (or any data structure)
#         3. Pass these values to one of the filters that will compute the necessary mathematical operations (mean, median, etc.)
#         4. Save the results at (i,j) in the ouput image.
#         5. return the output image
#         """
#
#         r, c = self.image.shape
#         imagenew = np.zeros((r, c))
#         size_1 = self.filter_size
#         size_1 = int(0.5 * size_1 - 0.5)
#         pad_image = np.zeros((r + (2 * size_1), c + (2 * size_1)))
#         pad_newimage = np.zeros((r + (2 * size_1), c + (2 * size_1)))
#         rx , ry = pad_image.shape
#         ker = []
#
#         for i in range(c):
#             for j in range(r):
#                 pad_image[i + size_1][ j + size_1] = self.image[i][ j]
#
#         for i in range(size_1 + 1, rx - size_1):
#             for j in range(size_1 + 1, rx - size_1):
#                 ker.clear()
#                 for x in range(-size_1, size_1 + 1):
#                     for y in range(-size_1, size_1 + 1):
#                         ker.append(pad_image[i + x][ j + y])
#                         if self.filter_name == 'alpha_trimmed':
#                             mask = self.filter(ker, self.alpha_d)
#                         elif self.filter_name == 'contra_harmonic':
#                             mask = self.filter(ker, self.order)
#                         else:
#                             mask = self.filter(ker)
#                         pad_newimage[i][j] = mask
#         for x1 in range(r):
#             for y1 in range(c):
#                 imagenew[x1][y1] = pad_newimage[x1 + size_1][y1 + size_1]
#         return imagenew


