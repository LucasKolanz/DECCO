class Param:
    '''
    Acts like a dictionary of dictionaries with default values for keys of
    the inner dictionaries.
    '''

    def __init__(self, defaults, methods, labels={}):
        self.adjusted = {m: {'method': m} for m in methods}
        for m in methods:
            self.adjusted[m].update(defaults)

        # Save some copy/paste for commonly used parameters
        self.set_param('todd', 'algorithm', 'wa')
        self.set_param('newton2', 'method', 'newton')
        self.set_param('newton2', 'initialize', 'ellipsoids')
        self.set_param('newtonr', 'method', 'newton')
        self.set_param('newtonr', 'initialize', 'random')
        self.set_param('bfgs', 'method', 'BFGS')
        self.set_param('lbfgs', 'method', 'L-BFGS')
        self._add_bfgs()

        # Save some copy/paste for commonly used labels
        self.labels = {'newton': 'Newton', 'todd': 'Coordinate Ascent',
                       'gradient': 'Gradient', 'newton2': 'Newton ellipsoids',
                       'cg': 'CG', 'truncated': 'Truncated Newton',
                       'newtonr': 'Newton random', 'bfgs': 'BFGS',
                       'lbfgs': 'L-BFGS'}
        self._add_bfgs_labels()
        for l in labels.keys():
            self.labels[l] = labels[l]


    def _add_bfgs(self):
        # More defaults for BFGS to save typing
        self.set_param('bfgs', 'method', 'BFGS')
        self.set_param('bfgs', 'bfgs_restart', 'update')
        self.set_param('bfgs', 'bfgs_method', 'BFGS')
        self.set_param('bfgsub', 'method', 'BFGS')
        self.set_param('bfgsub', 'bfgs_restart', 'update')
        self.set_param('bfgsub', 'bfgs_method', 'BFGS')
        self.set_param('bfgsrb', 'method', 'BFGS')
        self.set_param('bfgsrb', 'bfgs_restart', 'restart')
        self.set_param('bfgsrb', 'bfgs_method', 'BFGS')
        self.set_param('bfgsud', 'method', 'BFGS')
        self.set_param('bfgsud', 'bfgs_restart', 'update')
        self.set_param('bfgsud', 'bfgs_method', 'DFP')
        self.set_param('bfgsrd', 'method', 'BFGS')
        self.set_param('bfgsrd', 'bfgs_restart', 'restart')
        self.set_param('bfgsrd', 'bfgs_method', 'DFP')
        self.set_param('bfgsubt', 'method', 'BFGS')
        self.set_param('bfgsubt', 'bfgs_restart', 'update')
        self.set_param('bfgsubt', 'bfgs_method', 'BFGS')
        self.set_param('bfgsubt', 'search', 'trust')
        self.set_param('bfgsrbt', 'method', 'BFGS')
        self.set_param('bfgsrbt', 'bfgs_restart', 'restart')
        self.set_param('bfgsrbt', 'bfgs_method', 'BFGS')
        self.set_param('bfgsrbt', 'search', 'trust')


    def _add_bfgs_labels(self):
        # More defaults for BFGS to save typing
        labels = {'bfgsrb': 'BFGS-R-B', 'bfgsrd': 'BFGS-R-D',
                  'bfgsub': 'BFGS-U-B', 'bfgsud': 'BFGS-U-D'}
        self.labels.update(labels)


    def get_param(self, method, key):
        return self.adjusted[method][key]


    def get_label(self, method):
        try:
            label = self.labels[method]
        except KeyError:
            label = method
        return label


    def get_dict(self, method):
        return self.adjusted[method]


    def set_param(self, method, key, value):
        try:
            self.adjusted[method][key] = value
        except KeyError:
            pass
