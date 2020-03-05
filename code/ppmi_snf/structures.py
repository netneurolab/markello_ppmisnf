# -*- coding: utf-8 -*-
"""
Data objects and supporting functions for "easily" holding mixed project data
in HDF5 files
"""

import os
import os.path as op

import h5py
import numpy as np
import pandas as pd


def _recursive_save(h5file, obj, group):
    """
    Recursively saves `obj` to `h5file` in `group`

    Parameters
    ----------
    h5file : :obj:`h5py.File`
    obj : dict
    group : str
        Group in `h5file` in which to create datasets
    """

    grp = h5file.create_group(group) if group not in h5file else h5file[group]
    for key, item in obj.items():
        if isinstance(item, dict):
            _recursive_save(h5file, item, group=group + '/' + key)
        elif isinstance(item, np.ndarray):
            if key in grp:
                del grp[key]
            grp.create_dataset(key, item.shape, item.dtype)[...] = item
        elif item is not None:
            coerce = np.string_ if isinstance(item, str) else np.asarray
            item = coerce(item)
            if key in grp:
                del grp[key]
            grp.create_dataset(key, item.shape, item.dtype)[...] = item


def _recursive_load(h5file, group):
    """
    Recursively loads data from `h5file`

    Parameters
    ----------
    h5file : :obj:`h5py.File`
    group : str
        Group in `h5file` from which to load datasets

    Returns
    -------
    results : dict
        Dictionary containing loaded data
    """

    results = dict()
    for key, item in h5file[group].items():
        if isinstance(item, h5py.Dataset):
            results[key] = item[()]
            if hasattr(results[key], 'decode'):
                results[key] = results[key].decode()
        elif isinstance(item, h5py.Group):
            results[key] = _recursive_load(h5file, group=group + '/' + key)

    return results


def _get_keys(h5file, group='/'):
    """
    Returns list of all keys (datasets) in `h5file`

    Parameters
    ----------
    h5file : h5py.File
    group : str, optional
        Starting group in `h5file` from which to check for datasets.

    Returns
    -------
    keys : set
        Groups identified in `h5file`
    """

    keys = set()
    for key, item in h5file[group].items():
        not_pandas = h5file[group].attrs.get('df') is None
        if isinstance(item, h5py.Group) and not_pandas:
            keys.update(_get_keys(h5file, group=(group + key + '/')))
        else:
            keys.add(group + key)

    return keys


def _get_groups(h5file, group='/'):
    """
    Returns list of all groups in `h5file`
    """

    groups = set()
    for key, item in h5file[group].items():
        if isinstance(item, h5py.Group):
            groups.add(group[:-1])
            if h5file[group].attrs.get('df') is None:
                groups.add(group + key)
                groups.update(_get_groups(h5file, group=(group + key + '/')))
    return groups


class Frog():
    """
    A small helper class for storing mixed project data in an HDF5 file

    Named after Frog, my dog, who will eat just about anything
    """

    # map specifying valid data type(s) to store in group
    _map = {
        '/raw': (pd.DataFrame,),
        '/regressors': (pd.DataFrame,),
        '/processed': (pd.DataFrame,),
        '/snf': (np.ndarray, dict)
    }

    def __init__(self, filename):
        filename = op.abspath(filename)
        if not op.exists(filename):
            os.makedirs(op.dirname(filename), exist_ok=True)
            with h5py.File(filename, 'a'):
                pass
        elif not h5py.is_hdf5(filename):
            raise TypeError(f'Provided file {filename} exists but it is not a '
                            'valid HDF5 file.')

        self.filename = op.abspath(filename)
        self._mode = 'r'
        self._savers = {
            pd.DataFrame: self._save_dataframe,
            np.ndarray: self._save_array,
            dict: self._save_dict
        }
        self._loaders = {
            pd.DataFrame: self._load_dataframe,
            np.ndarray: self._load_array,
            dict: self._load_dict
        }

    def __enter__(self):
        self._fobj = h5py.File(self.filename, self._mode)
        return self._fobj

    def __exit__(self, type, value, tb):
        self._fobj.close()
        del self._fobj

    def __delitem__(self, key):
        with self._open('r+') as h5file:
            if key in h5file:
                del h5file[key]

    def __getitem__(self, key):
        return self.load(key)

    def __setitem__(self, key, value):
        return self.save(value, key, overwrite=True)

    def __contains__(self, key):
        return key in self.groups() + self.keys()

    def _ipython_key_completions_(self):
        return self.groups() + self.keys()

    def _open(self, mode='r'):
        """ Sets mode of object to open
        """
        self._mode = mode
        return self

    def _check_same(self, obj, group, compfunc):
        """ Checks whether `group` exists and contents are identical to `obj`

        Returns
        -------
        exists : bool
            Whether `group` exists
        same : bool
            Whether contents of `group` are identical to `obj`
        """
        if self._check_key(group):
            exists = same = True
            try:
                compfunc(obj, self.load(group))
            except AssertionError:
                same = False
        else:
            exists = same = False

        return exists, same

    def _save_dict(self, obj, group, overwrite=False):
        """ Saves dictionary `obj` to disk stored at `group`
        """
        compfunc = np.testing.assert_equal
        exists, same = self._check_same(obj, group, compfunc)

        if (not same and overwrite) or not exists:
            with self._open('r+') as h5file:
                _recursive_save(h5file, obj, group)

        return self

    def _load_dict(self, group):
        """ Returns dictionary stored at `group`
        """
        with self._open('r') as h5file:
            d = _recursive_load(h5file, group)
        return d

    def _save_array(self, arr, group, overwrite=False):
        """ Saves array `arr` to disk stored at `group`
        """
        compfunc = np.testing.assert_array_almost_equal
        exists, same = self._check_same(arr, group, compfunc)

        if (not same and overwrite) or not exists:
            group, key = op.split(group)
            with self._open('r+') as h5file:
                grp = h5file.get(group)
                if grp is None:
                    grp = h5file.create_group(group)
                grp.create_dataset(key, arr.shape, arr.dtype)[...] = arr

        return self

    def _load_array(self, group):
        """ Returns array stored at `group`
        """
        with self._open('r') as h5file:
            a = h5file[group][()]
        return a

    def _save_dataframe(self, df, group, overwrite=False):
        """ Saves dataframe `df` to disk stored at `group`
        """
        compfunc = pd.testing.assert_frame_equal
        exists, same = self._check_same(df, group, compfunc)

        if (not same and overwrite) or not exists:
            df.to_hdf(self.filename, group, mode='r+', format='table')
            with self._open('r+') as h5file:
                h5file[op.split(group)[0]].attrs['df'] = 1

        return self

    def _load_dataframe(self, group):
        """ Returns dataframe stored at `group`
        """
        return pd.read_hdf(self.filename, key=group)

    def _get_dtypes(self, group):
        """ Returns expected datatype for `group`
        """

        grp = key = group
        while key != '':
            for k, d in self.__class__._map.items():
                if grp == k:
                    return d
            grp, key = op.split(grp)

        raise ValueError(f'Invalid group {group}. Must be subset of '
                         f'{list(self.__class__._map.keys())}')

    def _check_key(self, key):
        """ Returns True if `key` exists; otherwise returns False
        """
        exists = True
        with self._open('r') as h5file:
            if h5file.get(key) is None:
                exists = False
        return exists

    def save(self, data, group, overwrite=False):
        """ Writes `data` to `group` with checks for datatype
        """
        if not group.startswith('/'):
            group = '/' + group
        dtypes = self.__class__._map.get('/' + group.split('/')[1])

        # determine expected datatype of `group`
        dtype = None
        for d in self._get_dtypes(group):
            if isinstance(data, d):
                dtype = d
        if dtype is None:
            raise TypeError(f'Cannot save data of type {type(data)} to '
                            f'{group}. Must be of type: {dtypes}.')

        return self._savers.get(dtype)(data, group, overwrite=overwrite)

    def load(self, group='/'):
        """ Loads and returns `group`
        """
        if group == '/':  # get everything in the file
            return {key: self.load(key) for key in self.__class_._map.keys()}

        if not group.startswith('/'):
            group = '/' + group
        if group.endswith('/'):
            group = group[:-1]

        if group in self.groups():      # it's a group
            grp, key = group, ''
        elif group in self.keys():      # it's a dataset
            grp, key = op.split(group)
        else:                           # it's invalid
            raise KeyError(f'{group} does not exist in file.')

        dtypes = self._get_dtypes(grp)
        with self._open('r') as h5file:
            isgrp = isinstance(h5file[group], h5py.Group)
            keys = sorted(list(_get_keys(h5file, grp + '/')))

        if dict in dtypes and isgrp:
            keys += [group]
            loader = self._loaders.get(dict)
        else:
            loader = self._loaders.get(dtypes[0])

        if group in keys:
            out = loader(group)
        else:
            out = {op.relpath(key, grp): loader(key) for key in keys}

        return out

    def keys(self):
        """ Lists keys of all datasets stored in file
        """
        k = []
        for grp in self.__class__._map.keys():
            if self._check_key(grp):
                with self._open('r') as h5file:
                    grpkeys = _get_keys(h5file, grp + '/')
                for key in grpkeys:
                    if self._check_key(key):
                        k.append(key)
        return sorted(list(k))

    def groups(self):
        """ Lists all groups stored in file
        """
        g = []
        for grp in self.__class__._map.keys():
            if self._check_key(grp):
                with self._open('r') as h5file:
                    g.extend(_get_groups(h5file, grp + '/'))
        return sorted(list(g))
