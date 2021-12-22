import matplotlib.pyplot as plt
import numpy
import numpy as np
import math as m
import cv2


def my_fun(x, *args):
    """
    my_fun generates the desired profile for the branch
    :param x: float with x-coordinate
    :return: value of function
    """
    br, r, delta = args[0], args[1], args[2]

    # Check the inputs
    assert r < 0.5, 'The radius has to be less than 0.5'
    assert 0 <= x <= 1, 'x has to be contained between 0 and 1'
    assert r < br, 'br has to be greater than r'

    if delta <= x <= 1-delta:
        # Scale x
        x = (x - delta)/(1 - 2*delta)
        # Evaluate the function
        if 0 <= x < br-r:
            return (1 - r) / (br - r) * x
        elif br-r <= x < br:
            return m.sqrt(r**2 - (x - br)**2) + (1 - r)
        elif br <= x < 1-br:
            return 1
        elif 1-br <= x < 1-br+r:\
            return m.sqrt(r**2 - (x - 1 + br)**2) + (1 - r)
        elif 1-br+r <= x <= 1:
            return (1 - r) / (br - r) * (1 - x)
    else:
        return 0


def evaluate_fourier(a, b, x):
    """
    evaluate_fourier gets the value of the fourier series for x, given the
    coefficients of the series.
    :param a: array of coefficients which multiply cos(x)
    :param b: array of coefficients which multiply sin(x)
    :param x: point where we want to evaluate the series
    :return: value of the sum
    """
    assert 0 <= x <= 1, "x variable not contained in [0, 1] \n"

    # Generate the sum for both arrays
    sum_a = sum((a[k] * m.cos(2 * m.pi * (k) * x) for k in range(len(a))))
    sum_b = sum((b[k] * m.sin(2 * m.pi * (k) * x) for k in range(len(b))))

    # Return the output
    return sum_a + sum_b - a[0] / 2


def generate_fourier(a, b, base_radius, branches, impact, n):
    image_pixel_size = 2048 * 2
    image_size = 0.2
    ext_radius = 0.07
    image = np.zeros((image_pixel_size, image_pixel_size), np.uint8)
    meters_per_pixel = image_size / image_pixel_size
    pixel_base_radius = m.floor(base_radius / meters_per_pixel)  # Get the base radius in pixels

    # Build the polygon
    points = []
    for i in range(branches):
        for k in range(n):
            modifier = 1 + impact * evaluate_fourier(a, b, k / n)
            points.append(compute_point_image_position(image_pixel_size, pixel_base_radius,
                                                       branches, i, k, n, modifier))
    points = np.floor(np.array(points))

    # Draw the shape
    cv2.fillPoly(image, np.int32([points]), 1, 255)
    fig = plt.figure(facecolor='w', figsize=(15, 15))

    fig.suptitle('Geometry drawing')
    plt.imshow(image, cmap='gray')


def generate_fourier_coefficients(n_coefs, period, fun, *args):
    """ generate_fourier_coefficients for given configuration """

    # Generate the x coordinates
    x = np.linspace(0, 1, 1000)

    # Preallocate results as numpy arrays
    a_s = np.empty(shape=(n_coefs, 1))
    b_s = np.empty(shape=(n_coefs, 1))

    # equations 13 onwards used from https://mathworld.wolfram.com/FourierSeries.html
    for i in range(0, n_coefs):
        sin_fourier = lambda x_: m.sin(2 * m.pi * i * x_ / period) * fun(x_, *args)
        cos_fourier = lambda x_: m.cos(2 * m.pi * i * x_ / period) * fun(x_, *args)
        sin_fourier = np.vectorize(sin_fourier)
        cos_fourier = np.vectorize(cos_fourier)
        a_s[i] = 2 / period * np.trapz(cos_fourier(x), x)
        b_s[i] = 2 / period * np.trapz(sin_fourier(x), x)

    # Return the results
    return a_s, b_s


def compute_point_image_position(image_pixel_size, pixel_base_radius, branches,
                                 i, k, n, modifier):
    return [image_pixel_size // 2 + pixel_base_radius * m.sin(2 * m.pi / branches * (i + k / n)) * modifier,
            image_pixel_size // 2 + pixel_base_radius * m.cos(2 * m.pi / branches * (i + k / n)) * modifier]


def generate_image(base_radius, branches, impact, n, br, r, delta):
    image_pixel_size = 2048*2
    image_size = 0.2
    ext_radius = 0.07
    image = np.zeros((image_pixel_size, image_pixel_size), np.uint8)
    meters_per_pixel = image_size/image_pixel_size
    pixel_base_radius = m.floor(base_radius / meters_per_pixel)  # Get the base radius in pixels

    # Build the polygon
    points = []
    for i in range(branches):
        for k in range(n):
            modifier = 1 + impact * my_fun(k / n, br, r, delta)
            points.append(compute_point_image_position(image_pixel_size, pixel_base_radius,
                                                       branches, i, k, n, modifier))
    points = np.floor(np.array(points))

    # Draw the shape
    cv2.circle(image, (image_pixel_size // 2, image_pixel_size // 2), int(ext_radius / meters_per_pixel), 255, 4)
    cv2.fillPoly(image, np.int32([points]), 255, 255)
    fig = plt.figure(facecolor='w', figsize=(15, 15))

    fig.suptitle('Geometry drawing')
    plt.imshow(image, cmap='gray')


if __name__ == '__main__':
    impact = 0.8
    br, r, delta = 0.3, 0.2, 0.17
    base_radius, branches, n = 0.032, 12, 50
    #generate_image(base_radius, branches, impact, n, br, r, delta)
    a, b = generate_fourier_coefficients(100, 1, my_fun, br, r, delta)

    x = np.linspace(0, 1, 100)
    y = []
    for i in range(len(x)):
        y.append(evaluate_fourier(a, b, x[i]))
    plt.plot(x, y)
    plt.grid()

    generate_fourier(a, b, base_radius, branches, impact, n)

    plt.show()


