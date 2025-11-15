import numpy as np
from scipy.signal import convolve2d

def create_flat_field(images):
    """
    Creates a CCD flat field from a set of spatially displaced data images
    using the Kuhn, Lin and Loranz method (Pub. Astron. Soc. Pacific, 103:1097 (1991)).

    Parameters
    ----------
    images : list of numpy.ndarray
        A list of 2D numpy arrays, where each array is a data image.
        Assumes all images have the same dimensions.

    Returns
    -------
    numpy.ndarray
        The calculated flat field image.
    """
    if not images:
        raise ValueError("The 'images' list cannot be empty.")

    # Convert images to float for calculations
    images = [img.astype(float) for img in images]

    # Get image dimensions
    ny, nx = images[0].shape
    num_images = len(images)

    # Step 1: Compute the median of all images to get an initial estimate of the flat
    # This step is common for most flat-fielding methods to get a rough idea of the scene
    median_image = np.median(np.array(images), axis=0)

    # Step 2: Compute the ratio of each image to the median image
    # These ratios will highlight the variations due to the flat field
    ratios = [img / median_image for img in images]

    # Step 3: Compute the mean of the ratios for each pixel
    # This mean represents the average non-uniformity across all images
    mean_ratios = np.mean(np.array(ratios), axis=0)

    # Step 4: Normalize the mean ratios by their average
    # This step ensures the flat field has an average value close to 1,
    # which is characteristic of a flat field (pixel values reflect sensitivity relative to the average)
    flat_field = mean_ratios / np.mean(mean_ratios)

    return flat_field

def apply_flat_field(data_image, flat_field):
    """
    Applies the calculated flat field to a data image.

    Parameters
    ----------
    data_image : numpy.ndarray
        The 2D data image to be corrected.
    flat_field : numpy.ndarray
        The calculated flat field.

    Returns
    -------
    numpy.ndarray
        The flat-field corrected image.
    """
    return data_image / flat_field

if __name__ == '__main__':
    # Example Usage:
    # Let's simulate some spatially displaced images with a non-uniform flat field.

    # Define a simulated true flat field (e.g., brighter in the center, darker at edges)
    ny_sim, nx_sim = 100, 100
    x = np.linspace(-1, 1, nx_sim)
    y = np.linspace(-1, 1, ny_sim)
    X, Y = np.meshgrid(x, y)
    true_flat = 1.0 - 0.3 * (X**2 + Y**2)
    true_flat = np.clip(true_flat, 0.5, 1.0) # Ensure values are positive

    print("Simulated True Flat Field:")
    print(true_flat[::10, ::10]) # Print a subset for brevity

    # Create a simulated "sky" background or object that shifts
    base_signal = 1000 * np.ones((ny_sim, nx_sim)) # A uniform background
    # Add a "hot" pixel or small feature for realism
    base_signal[ny_sim // 2, nx_sim // 2] += 500

    # Simulate multiple spatially displaced images
    num_sim_images = 5
    simulated_images = []
    for i in range(num_sim_images):
        # Simulate slight spatial displacement (e.g., 1-2 pixels)
        shift_x = int(np.random.randint(-2, 3))
        shift_y = int(np.random.randint(-2, 3))

        shifted_signal = np.roll(base_signal, shift_y, axis=0)
        shifted_signal = np.roll(shifted_signal, shift_x, axis=1)

        # Apply the true flat field and add some noise
        # The flat field is multiplicative
        noisy_image = shifted_signal * true_flat + np.random.normal(0, 5, (ny_sim, nx_sim))
        simulated_images.append(noisy_image)

    print(f"\nGenerated {num_sim_images} simulated images.")
    # Show a snippet of one simulated image
    print("Example of a simulated image (top-left corner):")
    print(simulated_images[0][:5, :5])

    # Calculate the flat field using the Kuhn, Lin and Loranz method
    calculated_flat = create_flat_field(simulated_images)

    print("\nCalculated Flat Field (top-left corner):")
    print(calculated_flat[:5, :5])

    print("\nDifference between True and Calculated Flat (top-left corner):")
    print((true_flat - calculated_flat)[:5, :5])

    # Visualize the results
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    im0 = axes[0].imshow(true_flat, cmap='hot', origin='lower')
    axes[0].set_title('Simulated True Flat Field')
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(calculated_flat, cmap='hot', origin='lower')
    axes[1].set_title('Calculated Flat Field')
    fig.colorbar(im1, ax=axes[1])

    im2 = axes[2].imshow(np.abs(true_flat - calculated_flat), cmap='viridis', origin='lower')
    axes[2].set_title('Absolute Difference (True vs. Calculated)')
    fig.colorbar(im2, ax=axes[2])

    plt.tight_layout()
    plt.show()

    # Apply the flat field to one of the simulated images
    original_image = simulated_images[0]
    corrected_image = apply_flat_field(original_image, calculated_flat)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    im3 = axes[0].imshow(original_image, cmap='gray', origin='lower')
    axes[0].set_title('Original Simulated Image (with flat field)')
    fig.colorbar(im3, ax=axes[0])

    im4 = axes[1].imshow(corrected_image, cmap='gray', origin='lower')
    axes[1].set_title('Flat-Field Corrected Image')
    fig.colorbar(im4, ax=axes[1])

    plt.tight_layout()
    plt.show()

    # You can also visually inspect the profiles
    plt.figure(figsize=(10, 5))
    plt.plot(true_flat[ny_sim // 2, :], label='True Flat Field (Middle Row)')
    plt.plot(calculated_flat[ny_sim // 2, :], label='Calculated Flat Field (Middle Row)')
    plt.title('Horizontal Profile of Flat Fields (Middle Row)')
    plt.xlabel('X-pixel')
    plt.ylabel('Relative Intensity')
    plt.legend()
    plt.grid(True)
    plt.show()