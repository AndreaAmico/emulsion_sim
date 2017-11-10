#python 3.5

import collections
import inspect
import lmfit
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import scipy.signal as signal


fit_plot_params = {
    'FILL_BETWEEN_COLOR' : [x/255 for x in (223, 180, 175)],
    'BACKGROUND_COLOR' : [x/255 for x in (254, 251, 243, 255)],
    'MARKER_COLOR' : (0.35,0.77,0.55),
    'BASE_COLOR' : 'brown',
    'LINEWIDTH' : 1.5,
    'PLOT_POINTS' : 100,
    'FIGSIZE' : (8,5),
    'PADDING' : 0.1,
    'MARKER' : 'o',
    'MARKEREDGEWIDTH_BIG' : 0.8,
    'MARKEREDGEWIDTH_SMALL' : 0.6,
    'MARKERSIZE_BIG' : 7,
    'MARKERSIZE_SMALL' : 4,
    'NUM_POINT_LIMIT' : 100,
    'GRID' : True,
    'GRID_LINEWIDTH' : 0.5,
    'GRID_STYLE' : ':'}


def confidence_interval(function, x, lmfit_output, confidence_probability=0.95):
    '''Returns the confidence interval of a fit
    
    Parameters
    ----------
    function : function
        Fitting model function(x, *params) = y
    x : list or numpy array
        list of coordinates where to evaluate the confidence interval on the variable y
    lmfit_output : lmfit.minimizer.MinimizerResult
        result of the fit
    confidence_probability : float    
    
    '''
    # best_fit_values = [lmfit_output.params[param].value for param in lmfit_output.params]
    not_fixed_params = [param for param in lmfit_output.params if lmfit_output.params[param].vary]
    
    # Compute the gradient of the fitting function with respect to the fitting parameters
    grad_array = np.zeros(len(not_fixed_params))

    sqrt_machine_precision = np.sqrt(np.finfo(type(1.0)).eps)
    for index, param in enumerate(not_fixed_params):
        start_param = lmfit_output.params.valuesdict().copy()
        end_param = lmfit_output.params.valuesdict().copy()
        start_param[param] += sqrt_machine_precision
        end_param[param] -= sqrt_machine_precision

        grad_array[index] = (function(x, **end_param)-function(x, **start_param))/(2*sqrt_machine_precision)

    # Student's t distribution
    t = np.abs(st.t.ppf(1-confidence_probability, lmfit_output.nfree))
    
    # Confidence interal estimation: grad(fun)_all_param.T * covariance_matrix * grad(fun)_all_param
    # reduced chi squared is considered to be 1
    return np.sqrt(np.dot(np.dot(lmfit_output.covar, grad_array), grad_array)) * t

def fit(x, y, function, plot=True, fig_ax=None, print_errors=True, **kwarg):
    '''Fit a function with lmfit module.
    Returns a <lmfit.minimizer.MinimizerResult> (see lmfit module).
    
    Parameters
    ----------
    
    x : list or numpy array
        x data points
    y : list or numpy array
        y data points
    function : function
        Fitting model function(x, *params) = y
        The default value of the function parameters changes the behaviour of the fit:
            - No default value: the parameter can vary in the fit and the intial guess is set to 1
            - Single number X: the parameter can vary in the fit and the initial guess is set to X
            - Sized X of len==2: x[1] must be the string "Fixed". X[0] is set as a fixed parameter
              and is not considered in the fit.
            - Sized X of len==3: the parameter can vary in the fit, the initial guess is set to X[0],
              The minimum and the maximum value the parameter can reach is set to X[1] and X[2] respectively.
    plot : bool
        If ``True`` print the best fit parameters, plot the dataset, the function evaluated with
        the best fit parameters and the confidence interval (default = 95%)
    fig_ax : list or tuple or None
        If ``None`` creates new matplotlib figure and axes
        fig_ax[0] = reference to an existing matplotlib figure
        fis_ax[1] = reference to an existing matplotlib axex 
    print_errors : bool
        If ``True`` print out the output message of the lmfit when it is not succeeded. It olso print out
        the reduced chi squared when it is bigger than 10.
    
    '''
    
    function_parameters = inspect.signature(function).parameters
    param_list = iter(function_parameters)
    
    lmfit_params = lmfit.Parameters()
    x_p = next(param_list)
    
    for param in param_list:
        current_param = function_parameters[param]

        # Setting defaults for initial parameters
        name = current_param.name 

        # The parameter is empty -> initial value is set to 1
        if current_param.default == inspect._empty:
            lmfit_params.add(name=name, value=1)
            continue

        # The parameter is an iterable ->
        # if len==2 the first value is a fixed value for the fitting function, the second one is not considered
        # if len==3 the first value is the guess parameter for the fit, the second and the third are the min and max boundaries
        if isinstance(current_param.default, collections.Sized):
            value = current_param.default[0]

            if len(current_param.default) == 2:
                error_msg = '''A fit parameter with two elements is considered fixed in the fit.
                The second element must be the string "Fixed"'''
                assert current_param.default[1].lower() == 'fixed', error_msg
                lmfit_params.add(name=name, value=value, vary=False)   
                continue
                
            if len(current_param.default) == 3:
                min_val = current_param.default[1]
                max_val = current_param.default[2]
                lmfit_params.add(name=name, value=value, min=min_val, max=max_val)    
                continue
            else:
                raise ValueError('Parameters can only be a number or an iterable of length 2 or 3')

        # The parameter is a number -> it is used as the initial gueess for the fit
        else:
            lmfit_params.add(name=name, value=current_param.default)
            continue

    def residual(lmfit_params, x, y):
        return (y - function(x, **{p.name:p.value for p in lmfit_params.values()}))

    # FIT
    out = lmfit.minimize(residual, lmfit_params, args=(x, y))
    
    if print_errors:
        if out.chisqr > 10:
            print('Warning: chi squared = {:.1f}'.format(out.chisqr))
        if not out.message == 'Fit succeeded.':
                print(out.message)
    
    if not plot:
        return out

    # PLOT
    if fig_ax:
        fig, ax = fig_ax
    else:
        fig, ax = plt.subplots(figsize=fit_plot_params['FIGSIZE'])

    #set the range of the plot
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    x_plot_range = (x_min - (x_max - x_min)*fit_plot_params['PADDING'], x_max + (x_max - x_min)*fit_plot_params['PADDING'])
    y_plot_range = (y_min - (y_max - y_min)*fit_plot_params['PADDING'], y_max + (y_max - y_min)*fit_plot_params['PADDING'])    
    ax.set_xlim(x_plot_range)
    ax.set_ylim(y_plot_range)
    
    x_plot = np.linspace(*x_plot_range, num=fit_plot_params['PLOT_POINTS'])

    if not 'Fit succeeded.' in out.message :
        ax.scatter(x, y)
        ax.plot(x_plot,  function(x_plot, **{p.name:p.value for p in lmfit_params.values()}), c="red")
        return out

    # Print the fitting results
    def round_sig(val, err, sig=2):
        if err == 0:
            return (0, 0)
        return (round(val, sig-int(np.floor(np.log10(abs(err))))-1),
                round(err, sig-int(np.floor(np.log10(abs(err))))-1))
    
    for param in out.params:
        val = out.params[param].value
        if out.params[param].vary:            
            err = out.params[param].stderr
            print('{name} = {} +/- {} ({percent:.2f}%)'.format(name=param, *round_sig(val, err, sig=2), percent=100*err/val))
        else:
            print('{} = {} (FIXED)'.format(param, val))
        
    # Plot the result fitting_function with the fitted parameters 
    y_plot =  function(x_plot, **{p.name:p.value for p in out.params.values()})
    ax.plot(x_plot,  y_plot, c=fit_plot_params['BASE_COLOR'], linewidth=fit_plot_params['LINEWIDTH'])
    
    # Plot 95% confidence interval
    y_error = np.array([confidence_interval(function, x, out, confidence_probability=0.95) for x in x_plot])

    ax.fill_between(x_plot, y_plot-y_error, y_plot+y_error, interpolate=True, color=fit_plot_params['FILL_BETWEEN_COLOR'])
    
    # Plot datapoints
    ax.plot(x, y, linestyle='None',marker=fit_plot_params['MARKER'], mfc=fit_plot_params['MARKER_COLOR'], mec=fit_plot_params['BASE_COLOR'],
             markersize = fit_plot_params['MARKERSIZE_BIG'] if len(x)<100 else fit_plot_params['MARKERSIZE_SMALL'],
             markeredgewidth = fit_plot_params['MARKEREDGEWIDTH_BIG'] if len(x)<fit_plot_params['NUM_POINT_LIMIT'] else fit_plot_params['MARKEREDGEWIDTH_SMALL'])

    fig.patch.set_facecolor(fit_plot_params['BACKGROUND_COLOR'])
    ax.set_facecolor(fit_plot_params['BACKGROUND_COLOR'])
    ax.grid(b=fit_plot_params['GRID'], linestyle=fit_plot_params['GRID_STYLE'], linewidth=fit_plot_params['GRID_LINEWIDTH'])
    
    return out



if __name__ == '__main__':
    
    def function(x, amp, x0=4, sigma=5, off=(0,'fixed')):
        return amp*np.exp(-(x-x0)**2/(2*sigma**2)) + off
    
    np.random.seed(42)
    x = np.random.random(30)*35
    y = function(x, amp=1, x0=5, sigma=7, off=-0.02) + (np.random.random(x.shape)-0.5)*0.2

    fit(x, y, function)
    plt.show()


def weighted_fit(x, y, w, function, plot=True, fig_ax=None, print_errors=True, **kwarg):
    '''Fit a function with lmfit module using weights.
    Returns a <lmfit.minimizer.MinimizerResult> (see lmfit module).
    
    Parameters
    ----------
    
    x : list or numpy array
        x data points
    y : list or numpy array
        y data points
    w : list or numpy array
        weights of the datapoints (standard error)
    function : function
        Fitting model function(x, *params) = y
        The default value of the function parameters changes the behaviour of the fit:
            - No default value: the parameter can vary in the fit and the intial guess is set to 1
            - Single number X: the parameter can vary in the fit and the initial guess is set to X
            - Sized X of len==2: x[1] must be the string "Fixed". X[0] is set as a fixed parameter
              and is not considered in the fit.
            - Sized X of len==3: the parameter can vary in the fit, the initial guess is set to X[0],
              The minimum and the maximum value the parameter can reach is set to X[1] and X[2] respectively.
    plot : bool
        If ``True`` print the best fit parameters, plot the dataset, the function evaluated with
        the best fit parameters and the confidence interval (default = 95%)
    fig_ax : list or tuple or None
        If ``None`` creates new matplotlib figure and axes
        fig_ax[0] = reference to an existing matplotlib figure
        fis_ax[1] = reference to an existing matplotlib axex 
    print_errors : bool
        If ``True`` print out the output message of the lmfit when it is not succeeded. It olso print out
        the reduced chi squared when it is bigger than 10.
    
    '''
    
    function_parameters = inspect.signature(function).parameters
    param_list = iter(function_parameters)
    
    lmfit_params = lmfit.Parameters()
    x_p = next(param_list)
    
    for param in param_list:
        current_param = function_parameters[param]

        # Setting defaults for initial parameters
        name = current_param.name 

        # The parameter is empty -> initial value is set to 1
        if current_param.default == inspect._empty:
            lmfit_params.add(name=name, value=1)
            continue

        # The parameter is an iterable ->
        # if len==2 the first value is a fixed value for the fitting function, the second one is not considered
        # if len==3 the first value is the guess parameter for the fit, the second and the third are the min and max boundaries
        if isinstance(current_param.default, collections.Sized):
            value = current_param.default[0]

            if len(current_param.default) == 2:
                error_msg = '''A fit parameter with two elements is considered fixed in the fit.
                The second element must be the string "Fixed"'''
                assert current_param.default[1].lower() == 'fixed', error_msg
                lmfit_params.add(name=name, value=value, vary=False)   
                continue
                
            if len(current_param.default) == 3:
                min_val = current_param.default[1]
                max_val = current_param.default[2]
                lmfit_params.add(name=name, value=value, min=min_val, max=max_val)    
                continue
            else:
                raise ValueError('Parameters can only be a number or an iterable of length 2 or 3')

        # The parameter is a number -> it is used as the initial gueess for the fit
        else:
            lmfit_params.add(name=name, value=current_param.default)
            continue

    def residual(lmfit_params, x, y, w):
        return (y - function(x, **{p.name:p.value for p in lmfit_params.values()}))/w

    # FIT
    out = lmfit.minimize(residual, lmfit_params, args=(x, y, w))
    
    if print_errors:
        if out.chisqr > 1000:
            print('Warning: chi squared = {:.1f}'.format(out.chisqr))
        if not out.message == 'Fit succeeded.':
                print(out.message)
    
    if not plot:
        return out

    # PLOT
    if fig_ax:
        fig, ax = fig_ax
    else:
        fig, ax = plt.subplots(figsize=fit_plot_params['FIGSIZE'])

    #set the range of the plot
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    x_plot_range = (x_min - (x_max - x_min)*fit_plot_params['PADDING'], x_max + (x_max - x_min)*fit_plot_params['PADDING'])
    y_plot_range = (y_min - (y_max - y_min)*fit_plot_params['PADDING'], y_max + (y_max - y_min)*fit_plot_params['PADDING'])    
    ax.set_xlim(x_plot_range)
    ax.set_ylim(y_plot_range)
    
    x_plot = np.linspace(*x_plot_range, num=fit_plot_params['PLOT_POINTS'])

    if not 'Fit succeeded.' in out.message :
        ax.errorbar(x, y, w, fmt='o')
        ax.plot(x_plot,  function(x_plot, **{p.name:p.value for p in lmfit_params.values()}), c="red")
        return out

    # Print the fitting results
    def round_sig(val, err, sig=2):
        if err == 0:
            return (0, 0)
        return (round(val, sig-int(np.floor(np.log10(abs(err))))-1),
                round(err, sig-int(np.floor(np.log10(abs(err))))-1))

    for param in out.params:
        val = out.params[param].value
        if out.params[param].vary:
            err = out.params[param].stderr
            print('{name} = {} +/- {} ({percent:.2f}%)'.format(name=param, *round_sig(val, err, sig=2), percent=100*err/val))
        else:
            print('{} = {} (FIXED)'.format(param, val))
        
    # Plot the result fitting_function with the fitted parameters 
    y_plot =  function(x_plot, **{p.name:p.value for p in out.params.values()})
    ax.plot(x_plot,  y_plot, c=fit_plot_params['BASE_COLOR'], linewidth=fit_plot_params['LINEWIDTH'])
    
    # Plot 95% confidence interval
    y_error = np.array([confidence_interval(function, x, out, confidence_probability=0.95) for x in x_plot])

    ax.fill_between(x_plot, y_plot-y_error, y_plot+y_error, interpolate=True, color=fit_plot_params['FILL_BETWEEN_COLOR'])
    
    # Plot datapoints
    ax.errorbar(x, y, w, linestyle='None', color=fit_plot_params['BASE_COLOR'] )
    ax.plot(x, y, linestyle='None',marker=fit_plot_params['MARKER'], mfc=fit_plot_params['MARKER_COLOR'], mec=fit_plot_params['BASE_COLOR'],
             markersize = fit_plot_params['MARKERSIZE_BIG'] if len(x)<100 else fit_plot_params['MARKERSIZE_SMALL'],
             markeredgewidth = fit_plot_params['MARKEREDGEWIDTH_BIG'] if len(x)<fit_plot_params['NUM_POINT_LIMIT'] else fit_plot_params['MARKEREDGEWIDTH_SMALL'])

    fig.patch.set_facecolor(fit_plot_params['BACKGROUND_COLOR'])
    ax.set_facecolor(fit_plot_params['BACKGROUND_COLOR'])
    ax.grid(b=fit_plot_params['GRID'], linestyle=fit_plot_params['GRID_STYLE'], linewidth=fit_plot_params['GRID_LINEWIDTH'])
    
    return out

##################################################################################################
######################################## GAUSSIAN FIT ############################################
##################################################################################################
def gauss_function(params, x):
    amp = params['amplitude'].value
    sigma = params['sigma'].value
    x0 = params['peak_position'].value
    off = params['offset'].value
    return off + amp * np.exp(-(x-x0)**2/(2*sigma**2))


def fit_gaussian(x_data, y_data, report=False, plot=False):
    """ NB: Requires numpy, lmfit,and matplotlib if plot=True
        function form: off + amp * np.exp(-(x-x0)**2/(2*sigma**2))
        variable name: ['amplitude', 'offset', 'peak_position', 'sigma']
        to extract parameters from the output use: out.params["amplitude"].value
    """

    data_array = np.array([x_data, y_data]).T
    x_sorted_data = data_array[np.argsort(data_array[:,0])]
    y_sorted_data = data_array[np.argsort(data_array[:,1])]

    min_y = y_sorted_data[0, :]
    max_y = y_sorted_data[-1, :]

    min_x = x_sorted_data[0, :]
    max_x = x_sorted_data[-1, :]


    if np.abs((min_x[1] + max_x[1])/2 - min_y[1]) < np.abs((min_x[1] + max_x[1])/2 - max_y[1]):
        amp_sign = 1
        peak_position = max_y[0]
        offset = min_y[1]
    else:
        amp_sign = -1
        peak_position = min_y[0]
        offset = max_y[1]

    amplitude = amp_sign*np.abs(max_y[1] - min_y[1])

    right_side = x_sorted_data[x_sorted_data[:, 0]>peak_position, :]
    left_side = x_sorted_data[x_sorted_data[:, 0]<peak_position, :]

    right_side[:, 1] = np.abs(right_side[:, 1] - offset - amplitude/2)
    right_fwhm = right_side[np.where(right_side[:, 1]==np.min(right_side[:, 1]))[0],0]

    left_side[:, 1] = np.abs(left_side[:, 1] - offset - amplitude/2)
    left_fwhm = left_side[np.where(left_side[:, 1]==np.min(left_side[:, 1]))[0],0]

    sigma = ((right_fwhm - left_fwhm) / (2*np.sqrt(2 * np.log(2))))[0]
    
    def gauss_function(params, x):
        amp = params['amplitude'].value
        sigma = params['sigma'].value
        x0 = params['peak_position'].value
        off = params['offset'].value
        return off + amp * np.exp(-(x-x0)**2/(2*sigma**2))

    def residual(params, x, data):
        model = gauss_function(params, x)
        return (data - model)

    params = lmfit.Parameters()
    params.add('amplitude', value=amplitude)
    params.add('offset', value=offset)
    params.add('peak_position', value=peak_position)
    params.add('sigma', value=sigma)

    out = lmfit.minimize(residual, params, args=(data_array[:, 0], data_array[:, 1]))

    if plot:
        x_plot = np.linspace(np.min(x_data), np.max(x_data), 100)
        plt.scatter(x_data, y_data)
        plt.plot(x_plot, gauss_function(out.params, x_plot), c="green")
        plt.xlabel("x_data")
        plt.ylabel("y_data")       

#         plt.savefig("out.svg")
        plt.show()

    if report: lmfit.report_fit(out, show_correl=False)
        
    return out


##################################################################################################
########################################### SINE FIT #############################################
##################################################################################################
def fit_sine(x_data, y_data, frequency_grid_size=100, max_min_estimation_size=5, report=True, plot=True):
    """ NB: Requires numpy, lmfit, scipy and matplotlib if plot=True
    Fit datapoint with a sine function of the form off+amp*sin(freq*x+phase)
    The initial guess of the frequency is obtained using the lombscargle algorithm
    Returns the output of lmfit.minimize
    """
    
    y_data_sorted = np.sort(y_data)
    y_data_min = np.mean(y_data_sorted[:max_min_estimation_size])
    y_data_max = np.mean(y_data_sorted[-max_min_estimation_size:])
    y_data_offset = (y_data_min + y_data_max)/2

    y_data_centred = y_data - y_data_offset

    total_time = np.max(x_data) - np.min(x_data)
    min_frequency = 1/total_time
    normval = x_data.shape[0]
    max_frequency = min_frequency*normval

    frequency_grid = np.linspace(min_frequency, max_frequency, frequency_grid_size)
        
    pgram = signal.lombscargle(x_data, y_data_centred, frequency_grid)

    lombscargle_spectrum = np.sqrt(4*(pgram/normval))
    maximum_indices = np.where(lombscargle_spectrum==np.max(lombscargle_spectrum))    
    central_frequency = frequency_grid[maximum_indices[0]]

    def sin_function(params, x):
        amp = params['amp'].value
        pshift = params['phase'].value
        freq = params['frequency'].value
        off = params['offset'].value
        
        return amp * np.sin(x * freq  + pshift) + off

    def residual(params, x, data):
        model = sin_function(params, x)
        return (data - model)

    params = lmfit.Parameters()
    params.add('amp', value=(y_data_max-y_data_min)/2.)
    params.add('offset', value=y_data_offset)
    params.add('phase', value=0.0, min=-np.pi/2., max=np.pi/2.)
    params.add('frequency', value=central_frequency, min=0)

    out = lmfit.minimize(residual, params, args=(x_data, y_data))
    
    if plot:
        fig, ax = plt.subplots(1,2, figsize=(12, 4))
        ax[0].scatter(frequency_grid, lombscargle_spectrum)
        ax[0].axvline(central_frequency, c="green") 
        ax[0].set_xlabel("x_data")
        ax[0].set_ylabel("Lomb Scargle spectrum")    

        x_plot = np.linspace(np.min(x_data), np.max(x_data), 100)
        ax[1].scatter(x_data, y_data)
        ax[1].plot(x_plot, sin_function(out.params, x_plot), c="green")
        ax[1].set_xlabel("x_data")
        ax[1].set_ylabel("y_data")       
        
        # fig.savefig("out.svg")
        plt.show()   

    if report: lmfit.report_fit(out, show_correl=False)

    return out



##################################################################################################
##################################### GAUSSIAN 2D FIT ############################################
##################################################################################################
def fit_gaussian_2d(data, angle=None, report=False, plot=False, initial_params=None):
    """ NB: Requires numpy, lmfit,and matplotlib if plot=Truie
    function form: off + amp * np.exp(-(xy[0]-x0)**2/(2*sx**2))*np.exp(-(xy[1]-y0)**2/(2*sy**2))
    variable name: ['amplitude', 'offset', 'x0', 'y0', 'sigma_x', 'sigma_y', ('angle')]
    to extract parameters from the output use: out.params["amplitude"].value        
    """
    if initial_params:
        x0 = initial_params["x0"].value
        y0 = initial_params["y0"].value
        sigma_x = initial_params["sigma_x"].value
        sigma_y = initial_params["sigma_y"].value
        offset = initial_params["offset"].value
        amplitude = initial_params["amplitude"].value
        
    else:
        out_x = fit_gaussian(np.arange(data.shape[1]), np.mean(data, 0), report=False, plot=False)
        out_y = fit_gaussian(np.arange(data.shape[0]), np.mean(data, 1), report=False, plot=False)

        x0 = out_x.params["peak_position"].value
        y0 = out_y.params["peak_position"].value

        sigma_x = out_x.params["sigma"].value
        sigma_y = out_y.params["sigma"].value

        offset = (out_x.params["offset"].value + out_y.params["offset"].value)/2
        amplitude = (out_x.params["amplitude"].value + out_y.params["amplitude"].value)/2

    if angle:
        def gauss_function_2d(params, xy):
            amp = params['amplitude'].value
            off = params['offset'].value
            sx = params['sigma_x'].value
            sy = params['sigma_y'].value
            x0 = params['x0'].value
            y0 = params['y0'].value
            t = np.deg2rad(params['angle'].value)
            
            x_exp = np.exp(-((xy[0]-x0)*np.cos(t)+(xy[1]-y0)*np.sin(t))**2/(2*sx**2))
            y_exp = np.exp(-((xy[1]-y0)*np.cos(t)+(xy[0]-x0)*np.sin(t))**2/(2*sy**2))
            return off + amp * x_exp * y_exp

    else:
        def gauss_function_2d(params, xy):
            amp = params['amplitude'].value
            off = params['offset'].value
            sx = params['sigma_x'].value
            sy = params['sigma_y'].value
            x0 = params['x0'].value
            y0 = params['y0'].value
            
            return off + amp * np.exp(-(xy[0]-x0)**2/(2*sx**2))*np.exp(-(xy[1]-y0)**2/(2*sy**2))

    def residual(params, xy, data):
        model = gauss_function_2d(params, xy)
        return (data - model)

    xy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))

    params = lmfit.Parameters()
    params.add('amplitude', value=amplitude)
    params.add('offset', value=offset)
    params.add('x0', value=x0)
    params.add('y0', value=y0)
    params.add('sigma_x', value=sigma_x)
    params.add('sigma_y', value=sigma_y)
    if angle: params.add('angle', value=angle)

    out = lmfit.minimize(residual, params, args=(xy, data))
    if plot:
        from mpl_toolkits import axes_grid1

        def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
            """Add a vertical color bar to an image plot."""
            divider = axes_grid1.make_axes_locatable(im.axes)
            width = axes_grid1.axes_size.AxesY(im.axes, aspect=1/aspect)
            pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
            current_ax = plt.gca()
            cax = divider.append_axes("right", size=width, pad=pad)
            plt.sca(current_ax)
            return im.axes.figure.colorbar(im, cax=cax, **kwargs)

        fig = plt.figure(figsize = (10,5))
        im = plt.imshow(data, cmap="viridis")
        add_colorbar(im)
        plt.contour(xy[0], xy[1], gauss_function_2d(out.params, xy), 10, colors='w', alpha=0.4)

        #plt.savefig("out.svg")
        plt.show()   
    if report: lmfit.report_fit(out, show_correl=False)
    return out

def gaussian_2d_volume(out):
    sx = out.params["sigma_x"].value
    sy = out.params["sigma_y"].value
    amp = out.params["amplitude"].value
    return 2*np.pi*sx*sy*amp
    


def progress_bar(current_value, max_value):
    progress = ((current_value+1)/max_value)*100
    if progress>98: progress=100
    print('\r[{0}{1}] {2:.1f}%'.format('#'*int(progress/2), ' '*(50-int(progress/2)), progress), end='')


def azimuthal_average(image, center, dr=1., sigmas = [1.,1.], weightMask = None, poisson=False): 
    y, x = np.indices(image.shape)
    r = np.hypot((x - center[0])*sigmas[1]/sigmas[0], (y - center[1]))

    nbins = int(np.round(r.max() / dr)+1)
    maxbin = nbins * dr
    bins = np.linspace(0, maxbin, nbins+1)
    bin_centers = (bins[1:]+bins[:-1])/2

    whichbin = np.digitize(r.flat, bins)
    if weightMask is None:
        weights = np.ones(image.shape)
    else:
        if not np.array_equal(image.shape, weightMask.shape): 
            raise ValueError('azimuthalAverage mask shape is not the same as the image shape: '+str(weightMask.shape)+ " vs. "+str(image.shape))
        else: weights = weightMask
    radial_profile = np.array([(image*weights).flat[whichbin==b].sum() / weights.flat[whichbin==b].sum() for b in np.arange(1,nbins)])
    if poisson:
        radial_var = np.array([np.sqrt(np.var((image*weights).flat[whichbin==b])) for b in np.arange(1,nbins)])
        return bin_centers[:-1], radial_profile, radial_var
    return bin_centers[:-1], radial_profile

def azimuthal_sum(image, center, dr=1., sigmas = [1.,1.], weightMask = None, poisson=False): 
    y, x = np.indices(image.shape)
    r = np.hypot((x - center[0])*sigmas[1]/sigmas[0], (y - center[1]))

    nbins = int(np.round(r.max() / dr)+1)
    maxbin = nbins * dr
    bins = np.linspace(0, maxbin, nbins+1)
    bin_centers = (bins[1:]+bins[:-1])/2

    whichbin = np.digitize(r.flat, bins)
    if weightMask is None:
        weights = np.ones(image.shape)
    else:
        if not np.array_equal(image.shape, weightMask.shape): 
            raise ValueError('azimuthalAverage mask shape is not the same as the image shape: '+str(weightMask.shape)+ " vs. "+str(image.shape))
        else: weights = weightMask
    radial_profile = np.array([(image*weights).flat[whichbin==b].sum() for b in np.arange(1,nbins)])
    if poisson:
        radial_var = np.array([np.sqrt(np.var((image*weights).flat[whichbin==b])) for b in np.arange(1,nbins)])
        return bin_centers[:-1], radial_profile, radial_var
    return bin_centers[:-1], radial_profile


def bin_data(data_x, data_y, bins):
    ''' Bin data along the x dirextion using bins as separator.
        It returns the data binned along the x an y direction plus the x and y error as the standard deviation
        of the mean on the datapoints inside the bin
        returns x, y, x_err, y_err as numpy array
    '''
    data_x = np.array(data_x)
    data_y = np.array(data_y)
    digitized = np.digitize(data_x, bins)
    x = np.array([data_x[digitized == i].mean() for i in range(1, len(bins)) if i in digitized])
    x_err = np.array([data_x[digitized == i].std() for i in range(1, len(bins)) if i in digitized])
    y = [data_y[digitized == i].mean() for i in range(1, len(bins)) if i in digitized]
    y_err = np.array([data_y[digitized == i].std()/np.sqrt(data_y[digitized == i].shape[0]) for i in range(1, len(bins)) if i in digitized])
    return x, np.array(y), x_err, y_err

def bin_data_err(data_x, data_y, err_y,  bins):
    ''' Bin data along the x dirextion using bins as separator.
        It returns the data binned along the x an y direction plus the x and y error as the standard deviation
        of the mean on the datapoints inside the bin. The y error obtained this way is square summed with the individual
        y_errors of the dataset.
    '''
    data_x = np.array(data_x)
    data_y = np.array(data_y)
    digitized = np.digitize(data_x, bins)
    x = np.array([data_x[digitized == i].mean() for i in range(1, len(bins)) if i in digitized])
    x_err = np.array([data_x[digitized == i].std() for i in range(1, len(bins)) if i in digitized])
    y = [data_y[digitized == i].mean() for i in range(1, len(bins)) if i in digitized]
    y_err = np.array([data_y[digitized == i].std()/np.sqrt(data_y[digitized == i].shape[0]) for i in range(1, len(bins)) if i in digitized])
    y_sum_var = np.array([np.sum(err_y[digitized == i]**2) for i in range(1, len(bins)) if i in digitized])
    y_err_tot = np.sqrt(y_err**2 + y_sum_var)
    return x, np.array(y), y_err_tot


def bin_image(image, binsize_x=1, binsize_y=None, aggregation_function=np.sum):
    sy, sx = image.shape
    if not binsize_y:
        binsize_y = binsize_x
        
    y_bins, x_bins = sy // binsize_y, sx // binsize_x
    crop = np.ogrid[(sy % binsize_y)//2: sy-((sy % binsize_y)//2+(sy%binsize_y)%2),
                   (sx % binsize_x)//2: sx-((sx % binsize_x)//2+(sx%binsize_x)%2)]
    cropped = image[crop]
    x_agg = aggregation_function(cropped.reshape(cropped.shape[0], x_bins, binsize_x), axis=2)
    return aggregation_function(x_agg.reshape(y_bins, binsize_y, x_agg.shape[1]), axis=1)


def bin_images(image, binsize_x=1, binsize_y=None, aggregation_function=np.sum):
    num_images, sy, sx = image.shape
    if not binsize_y:
        binsize_y = binsize_x
        
    y_bins, x_bins = sy // binsize_y, sx // binsize_x
    crop = np.ogrid[0:num_images,
                   (sy % binsize_y)//2: sy-((sy % binsize_y)//2+(sy%binsize_y)%2),
                   (sx % binsize_x)//2: sx-((sx % binsize_x)//2+(sx%binsize_x)%2)]
    cropped = image[crop]
    x_agg = aggregation_function(cropped.reshape(num_images, cropped.shape[1], x_bins, binsize_x), axis=3)
    return aggregation_function(x_agg.reshape(num_images,y_bins, binsize_y, x_agg.shape[2]), axis=2)

def play_bell():
    import winsound
    duration = 200  # millisecond
    freq = 440  # Hz
    for i in range(5):
        winsound.Beep(int(freq*(i/2+1)), duration)