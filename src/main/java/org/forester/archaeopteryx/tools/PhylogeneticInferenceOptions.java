// $Id:
// $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx.tools;

import java.io.File;

import org.forester.archaeopteryx.Configuration;
import org.forester.evoinference.distance.PairwiseDistanceCalculator.PWD_DISTANCE_METHOD;
import org.forester.msa.Mafft;

public final class PhylogeneticInferenceOptions {

    private static final int                 BOOTSTRAP_RESAMPLES_DEFAULT                  = 100;
    private static final PWD_DISTANCE_METHOD PWD_DISTANCE_METHOD_DEFAULT                  = PWD_DISTANCE_METHOD.KIMURA_DISTANCE;
    public static final long                 RANDOM_NUMBER_SEED_DEFAULT                   = 42L;
    private static final boolean             PERFORM_BOOTSTRAP_RESAMPLING_DEFAULT         = false;
    private static final double              msa_processing_max_allowed_gap_ratio_default = 0.5;
    private static final int                 msa_processing_min_allowed_length_default    = 50;
    private int                              _bootstrap_samples;
    private PWD_DISTANCE_METHOD              _pwd_distance_method;
    private long                             _random_number_generator_seed;
    private boolean                          _perform_bootstrap_resampling;
    private String                           _intermediate_files_base;
    private String                           _msa_prg_parameters;
    private boolean                          _execute_msa_processing;
    private boolean                          _msa_processing_remove_all_gap_columns;
    private double                           _msa_processing_max_allowed_gap_ratio;
    private int                              _msa_processing_min_allowed_length;
    private boolean                          _save_pwd_file;
    private boolean                          _save_processed_msa;
    private boolean                          _save_original_msa;
    private File                             _pwd_outfile;
    private File                             _processed_msa_outfile;
    private File                             _original_msa_outfile;

    public synchronized String getMsaPrgParameters() {
        return _msa_prg_parameters;
    }

    public synchronized void setMsaPrgParameters( final String msa_prg_parameters ) {
        _msa_prg_parameters = new String( msa_prg_parameters );
    }

    public synchronized String getIntermediateFilesBase() {
        return _intermediate_files_base;
    }

    public synchronized String getMsaPrg() {
        return "MAFFT";
    }

    public synchronized void setIntermediateFilesBase( final String intermediate_files_base ) {
        _intermediate_files_base = new String( intermediate_files_base );
    }

    public PhylogeneticInferenceOptions() {
        init();
    }

    // Deep copy.
    public synchronized PhylogeneticInferenceOptions copy() {
        final PhylogeneticInferenceOptions o = new PhylogeneticInferenceOptions();
        o._bootstrap_samples = _bootstrap_samples;
        o._pwd_distance_method = _pwd_distance_method;
        o._random_number_generator_seed = _random_number_generator_seed;
        o._perform_bootstrap_resampling = _perform_bootstrap_resampling;
        o._intermediate_files_base = new String( _intermediate_files_base );
        o._msa_prg_parameters = new String( _msa_prg_parameters );
        o._msa_processing_max_allowed_gap_ratio = _msa_processing_max_allowed_gap_ratio;
        o._msa_processing_min_allowed_length = _msa_processing_min_allowed_length;
        o._execute_msa_processing = _execute_msa_processing;
        o._msa_processing_remove_all_gap_columns = _msa_processing_remove_all_gap_columns;
        o._save_pwd_file = _save_pwd_file;
        o._save_processed_msa = _save_processed_msa;
        o._save_original_msa = _save_original_msa;
        if ( _pwd_outfile != null ) {
            o._pwd_outfile = new File( _pwd_outfile.toString() );
        }
        if ( _processed_msa_outfile != null ) {
            o._processed_msa_outfile = new File( _processed_msa_outfile.toString() );
        }
        if ( _original_msa_outfile != null ) {
            o._original_msa_outfile = new File( _original_msa_outfile.toString() );
        }
        return o;
    }

    private synchronized void init() {
        _bootstrap_samples = BOOTSTRAP_RESAMPLES_DEFAULT;
        _pwd_distance_method = PWD_DISTANCE_METHOD_DEFAULT;
        _random_number_generator_seed = RANDOM_NUMBER_SEED_DEFAULT;
        _perform_bootstrap_resampling = PERFORM_BOOTSTRAP_RESAMPLING_DEFAULT;
        _intermediate_files_base = "";
        _msa_prg_parameters = Mafft.getDefaultParameters();
        _msa_processing_max_allowed_gap_ratio = msa_processing_max_allowed_gap_ratio_default;
        _msa_processing_min_allowed_length = msa_processing_min_allowed_length_default;
        _execute_msa_processing = false;
        _msa_processing_remove_all_gap_columns = false;
        _save_pwd_file = false;
        _save_processed_msa = false;
        _save_original_msa = false;
        _pwd_outfile = null;
        _processed_msa_outfile = null;
        _original_msa_outfile = null;
    }

    public synchronized void setBootstrapSamples( final int bootstrap_samples ) {
        _bootstrap_samples = bootstrap_samples;
    }

    public synchronized int getBootstrapSamples() {
        return _bootstrap_samples;
    }

    public synchronized void setPwdDistanceMethod( final PWD_DISTANCE_METHOD pwd_distance_method ) {
        _pwd_distance_method = pwd_distance_method;
    }

    public synchronized PWD_DISTANCE_METHOD getPwdDistanceMethod() {
        return _pwd_distance_method;
    }

    public synchronized void setRandomNumberGeneratorSeed( final long random_number_generator_seed ) {
        _random_number_generator_seed = random_number_generator_seed;
    }

    public synchronized long getRandomNumberGeneratorSeed() {
        return _random_number_generator_seed;
    }

    public synchronized void setPerformBootstrapResampling( final boolean perform_bootstrap_resampling ) {
        _perform_bootstrap_resampling = perform_bootstrap_resampling;
    }

    public synchronized boolean isPerformBootstrapResampling() {
        return _perform_bootstrap_resampling;
    }

    public static PhylogeneticInferenceOptions createInstance( final Configuration configuration ) {
        final PhylogeneticInferenceOptions o = new PhylogeneticInferenceOptions();
        if ( configuration.getDefaultBootstrapSamples() >= 0 ) {
            o.setBootstrapSamples( configuration.getDefaultBootstrapSamples() );
        }
        return o;
    }

    public File getTempDir() {
        //TODO
        return new File( "/Users/zma/Desktop/tmp/" );
    }

    public void setMsaProcessingMaxAllowedGapRatio( final double msa_processing_max_allowed_gap_ratio ) {
        _msa_processing_max_allowed_gap_ratio = msa_processing_max_allowed_gap_ratio;
    }

    public double getMsaProcessingMaxAllowedGapRatio() {
        return _msa_processing_max_allowed_gap_ratio;
    }

    public void setMsaProcessingMinAllowedLength( final int msa_processing_min_allowed_length ) {
        _msa_processing_min_allowed_length = msa_processing_min_allowed_length;
    }

    public int getMsaProcessingMinAllowedLength() {
        return _msa_processing_min_allowed_length;
    }

    boolean isExecuteMsaProcessing() {
        return _execute_msa_processing;
    }

    void setExecuteMsaProcessing( final boolean execute_msa_processing ) {
        _execute_msa_processing = execute_msa_processing;
    }

    boolean isMsaProcessingRemoveAllGapColumns() {
        return _msa_processing_remove_all_gap_columns;
    }

    void setMsaProcessingRemoveAllGapColumns( final boolean msa_processing_remove_all_gap_columns ) {
        _msa_processing_remove_all_gap_columns = msa_processing_remove_all_gap_columns;
    }

    boolean isSavePwdFile() {
        return _save_pwd_file;
    }

    void setSavePwdFile( final boolean save_pwd_file ) {
        _save_pwd_file = save_pwd_file;
    }

    boolean isSaveProcessedMsa() {
        return _save_processed_msa;
    }

    void setSaveProcessedMsa( final boolean save_processed_msa ) {
        _save_processed_msa = save_processed_msa;
    }

    boolean isSaveOriginalMsa() {
        return _save_original_msa;
    }

    void setSaveOriginalMsa( final boolean save_original_msa ) {
        _save_original_msa = save_original_msa;
    }

    File getPwdOutfile() {
        return _pwd_outfile;
    }

    void setPwdOutfile( final File pwd_outfile ) {
        _pwd_outfile = pwd_outfile;
    }

    File getProcesseMsaOutfile() {
        return _processed_msa_outfile;
    }

    void setProcesseMsaOutfile( final File processed_msa_outfile ) {
        _processed_msa_outfile = processed_msa_outfile;
    }

    File getOriginalMsaOutfile() {
        return _original_msa_outfile;
    }

    void setOriginalMsaOutfile( final File original_msa_outfile ) {
        _original_msa_outfile = original_msa_outfile;
    }
}
