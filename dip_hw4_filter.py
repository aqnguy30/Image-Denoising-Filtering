"""dip_hw4_filter.py: Starter file to run howework 4"""

#Example Usage: ./dip_hw4_filter -i image
#Example Usage: python dip_hw4_filter.py -i image


__author__  = "Pranav Mantini"
__email__ = "pmantini@uh.edu"
__version__ = "1.0.0"

import cv2
import sys
from Denoise.Filtering import Filtering
from datetime import datetime
import numpy as np


def display_image(window_name, image):
    """A function to display image"""
    cv2.namedWindow(window_name)
    cv2.imshow(window_name, image)
    cv2.waitKey(0)

def get_bipolar_noise( image, noise_proba, noise_probb):
    """ Adds bipolar noise to the image
       takes as input:
       image: the input image
       noise_proba: probability of a pixel to be noisy
       returns a noisy image"""
    noisy_image = image.copy()
    rows, cols = image.shape

    for i in range(rows):
        for j in range(cols):
            n = np.random.random()
            if n < 0.5:
                n = np.random.random()
                if n < noise_proba:
                    noisy_image[i][j] = 0
                else:
                    noisy_image[i][j] = image[i][j]
            else:
                n = np.random.random()
                if n < noise_probb:
                    noisy_image[i][j] = 255
                else:
                    noisy_image[i][j] = image[i][j]

    return noisy_image

def get_gaussian_noise(image, mean, var):
    """ Adds gaussian noise to the image
               takes as input:
               image: the input image
               mean: gaussian distribution mean
               var: gaussian distribution variance
               returns a noisy image"""
    x, y = image.shape[0], image.shape[1]
    sigma = var ** 0.5
    gaussian = np.random.normal(mean, sigma, (x, y))
    noisy_image = image + gaussian
    return noisy_image


def main():
    """ The main funtion that parses input arguments, calls the approrpiate
     fitlering method and writes the output image"""

    #Parse input arguments
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument("-i", "--image", dest="image",
                        help="specify the name of the image", metavar="IMAGE")
    parser.add_argument("-n", "--noise", dest="noise",
                        help="specify type of the noise to be added (gaussian, bipolar)", metavar="NOISE")
    parser.add_argument("-m", "--mask", dest="mask",
                        help="specify name of the mask (median, max, min, alpha_trimmed, arithmetic_mean, geometric_mean, contra_harmonic)", metavar="MASK")
    parser.add_argument("-s", "--mask_size", dest="mask_size",
                        help="specify the size of the filter", metavar="MASK SIZE")
    parser.add_argument("-p", "--alpha_d", dest="alpha_d",
                        help="specify the trimming parameter for alpha trimmed filter", metavar="P")
    parser.add_argument("-o", "--order", dest="order",
                        help="specify the order parameter for contra harmonic filter", metavar="O")
    parser.add_argument("-npa", "--noise_proba", dest="noise_proba",
                        help="specify the probability of pepper (a) noise", metavar="NRA")
    parser.add_argument("-npb", "--noise_probb", dest="noise_probb",
                        help="specify the probability of salt (b) noise", metavar="NRB")
    parser.add_argument("-mean", "--mean", dest="mean",
                        help="specify the mean parameter for th gaussian noise", metavar="MEAN")
    parser.add_argument("-v", "--var", dest="var",
                        help="specify the variance parameter for th gaussian noise", metavar="VAR")

    args = parser.parse_args()

    #Load image
    if args.image is None:
        print("Please specify the name of image")
        print("use the -h option to see usage information")
        sys.exit(2)
    else:
        image_name = args.image.split(".")[0]
        input_image = cv2.imread(args.image, 0)
        rows, cols = input_image.shape

    if args.noise is None:
        print("Noise not specified using default (gaussian)")
        print("use the -h option to see usage information")
        noise = 'gaussian'
    else:
        noise = args.noise

    #Check resize scale parametes
    if args.mask is None:
        print("Mask not specified using default (median)")
        print("use the -h option to see usage information")
        mask = 'median'
    elif args.mask not in ['median', 'max', 'min', 'alpha_trimmed','arithmetic_mean','geometric_mean','contra_harmonic']:
        print("Unknown mask, using default (median)")
        print("use the -h option to see usage information")
        mask = 'median'
    else:
        mask = args.mask

    if args.mask_size is None:
        print("Mask size not specified using default (3)")
        print("use the -h option to see usage information")
        mask_size = 3
    else:
        mask_size = int(args.mask_size)

    if mask == 'alpha_trimmed':
        if args.alpha_d is None:
            print("trimming parameter not specified, using default (2)")
            print("use the -h option to see usage information")
            alpha_d = 2
        else:
            alpha_d = float(args.alpha_d)
    elif mask == 'contra_harmonic':
        if args.order is None:
            print("order parameter not specified, using default 1.5")
            print("use the -h option to see usage information")
            order = 1.5
        else:
            order = float(args.order)

    noisy_image = input_image.copy()

    if noise == 'bipolar':
        if args.noise_proba is None:
            print("amount of pepper noise not specified, using default (0.01)")
            print("use the -h option to see usage information")
            noise_proba = 0.01
        else:
            noise_proba = float(args.noise_proba)


        if args.noise_probb is None:
            print("amount of salt noise not specified, using default (0.01)")
            print("use the -h option to see usage information")
            noise_probb = 0.01
        else:
            noise_probb = float(args.noise_probb)

        noisy_image = get_bipolar_noise(input_image, noise_proba, noise_probb)

    if noise == 'gaussian':
        if args.mean is None:
            print("the mean for gaussian noise is not specified, using default mean=0")
            print("use the -h option to see usage information")
            mean = 20
        else:
            mean = float(args.mean)
        if args.var is None:
            print("the var for gaussian noise is not specified, using default var=0.1")
            print("use the -h option to see usage information")
            var = 50
        else:
            var = float(args.var)
        noisy_image = get_gaussian_noise(input_image, mean, var)

    if mask == 'alpha_trimmed':
        Filter_obj = Filtering(noisy_image, mask, mask_size, alpha_d = alpha_d)
        output = Filter_obj.filtering()
    elif mask == 'contra_harmonic':
        Filter_obj = Filtering(noisy_image, mask, mask_size, order=order)
        output = Filter_obj.filtering()
    else:
        Filter_obj = Filtering(noisy_image, mask, mask_size)
        output = Filter_obj.filtering()

    #Write output file
    output_dir = 'output/'

    output_image_name = output_dir+image_name+"_"+noise+"_noise_"+datetime.now().strftime("%m%d-%H%M%S")+".jpg"
    cv2.imwrite(output_image_name, noisy_image)
    output_image_name = output_dir + image_name+"_denoised_" + mask + datetime.now().strftime("%m%d-%H%M%S") + ".jpg"
    cv2.imwrite(output_image_name, output)



if __name__ == "__main__":
    main()







