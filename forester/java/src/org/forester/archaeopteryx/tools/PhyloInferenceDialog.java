// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.border.Border;
import javax.swing.border.LineBorder;

import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.evoinference.distance.PairwiseDistanceCalculator.PWD_DISTANCE_METHOD;
import org.forester.sequence.Sequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

public class PhyloInferenceDialog extends JDialog implements ActionListener {

    private static final long                  serialVersionUID = 8337543508238133614L;
    private final JPanel                       _pnl;
    private final JButton                      _launch_btn;
    private final JButton                      _cancel_btn;
    private final JFormattedTextField          _bootstrap_tf;
    private final JCheckBox                    _bootstrap_cb;
    private final PhylogeneticInferenceOptions _opts;
    private JTextField                         _input_msa_file_tf;
    private JButton                            _select_input_msa_btn;
    private final MainFrameApplication         _parent_frame;
    private JTextField                         _msa_length_tf;
    private JTextField                         _msa_size_tf;
    private JTextField                         _msa_type_tf;
    private final JRadioButton                 _distance_calc_kimura_rb;
    private final JRadioButton                 _distance_calc_poisson_rb;
    private final JRadioButton                 _distance_calc_fract_dissimilarity_rb;
    private int                                _value           = JOptionPane.CANCEL_OPTION;
    private JTextField                         _input_seqs_tf;
    private JButton                            _select_input_seqs_btn;
    private JTextField                         _input_seqs_number_tf;
    private JTextField                         _input_seqs_median_length_tf;
    private JTextField                         _input_seqs_min_length_tf;
    private JTextField                         _input_seqs_max_length_tf;
    private JTextField                         _input_seqs_type_tf;
    private JTextField                         _mafft_paramenters_tf;
    private JTextField                         _clustalo_paramenters_tf;
    private JTextField                         _msa_processing_max_allowed_gap_ratio_tf;
    private JTextField                         _msa_processing_min_allowed_length_tf;
    private JTextField                         _random_seed_tf;
    private JCheckBox                          _execute_msa_processing_cb;
    private JCheckBox                          _msa_processing_remove_all_gap_columns_cb;
    private JCheckBox                          _mafft_cb;
    private JCheckBox                          _clustalo_cb;
    private JCheckBox                          _save_pwd_file_cb;
    private JCheckBox                          _save_processed_msa_cb;
    private JCheckBox                          _save_original_msa_cb;
    private JTextField                         _pwd_outfile_tf;
    private JTextField                         _processed_msa_outfile_tf;
    private JTextField                         _original_msa_outfile_tf;

    public PhyloInferenceDialog( final MainFrameApplication frame,
                                 final PhylogeneticInferenceOptions options,
                                 final boolean from_unaligned_seqs ) {
        super( frame, true );
        setVisible( false );
        _parent_frame = frame;
        _opts = options;
        _pnl = new JPanel();
        getContentPane().add( _pnl );
        final BoxLayout box_layout = new BoxLayout( _pnl, BoxLayout.PAGE_AXIS );
        _pnl.setLayout( box_layout );
        if ( from_unaligned_seqs ) {
            setTitle( "Phylogenetic Inference (including multiple sequence alignment)" );
            final JPanel inputfile_pnl_1 = new JPanel();
            final JPanel inputfile_pnl_2 = new JPanel();
            final JPanel inputfile_pnl_3 = new JPanel();
            final JPanel inputfile_pnl_4 = new JPanel();
            inputfile_pnl_1.setLayout( new FlowLayout() );
            inputfile_pnl_2.setLayout( new FlowLayout() );
            inputfile_pnl_3.setLayout( new FlowLayout() );
            inputfile_pnl_4.setLayout( new FlowLayout() );
            inputfile_pnl_1.add( new JLabel( "Input Sequence File:" ) );
            inputfile_pnl_1.add( _input_seqs_tf = new JTextField() );
            inputfile_pnl_1.add( _select_input_seqs_btn = new JButton( "Select Input File" ) );
            inputfile_pnl_2.add( new JLabel( "Sequences: " ) );
            inputfile_pnl_2.add( new JLabel( "Number of Sequences:" ) );
            inputfile_pnl_2.add( _input_seqs_number_tf = new JTextField() );
            inputfile_pnl_2.add( new JLabel( "Length: median:" ) );
            inputfile_pnl_2.add( _input_seqs_median_length_tf = new JTextField() );
            inputfile_pnl_2.add( new JLabel( "min:" ) );
            inputfile_pnl_2.add( _input_seqs_min_length_tf = new JTextField() );
            inputfile_pnl_2.add( new JLabel( "max:" ) );
            inputfile_pnl_2.add( _input_seqs_max_length_tf = new JTextField() );
            inputfile_pnl_2.add( new JLabel( "Type:" ) );
            inputfile_pnl_2.add( _input_seqs_type_tf = new JTextField() );
            inputfile_pnl_3.add( _mafft_cb = new JCheckBox( "MAFFT" ) );
            inputfile_pnl_3.add( new JLabel( "Parameters: " ) );
            inputfile_pnl_3.add( _mafft_paramenters_tf = new JTextField() );
            inputfile_pnl_4.add( _clustalo_cb = new JCheckBox( "ClustalO" ) );
            inputfile_pnl_4.add( new JLabel( "Parameters: " ) );
            inputfile_pnl_4.add( _clustalo_paramenters_tf = new JTextField() );
            _input_seqs_median_length_tf.setColumns( 4 );
            _input_seqs_min_length_tf.setColumns( 4 );
            _input_seqs_max_length_tf.setColumns( 4 );
            _input_seqs_number_tf.setColumns( 4 );
            _input_seqs_type_tf.setColumns( 2 );
            _input_seqs_tf.setColumns( 20 );
            _input_seqs_tf.setEditable( false );
            _input_seqs_median_length_tf.setEditable( false );
            _input_seqs_min_length_tf.setEditable( false );
            _input_seqs_max_length_tf.setEditable( false );
            _input_seqs_number_tf.setEditable( false );
            _input_seqs_type_tf.setEditable( false );
            _mafft_paramenters_tf.setColumns( 26 );
            _mafft_paramenters_tf.setText( "--maxiterate 1000 --localpair" );
            _clustalo_paramenters_tf.setColumns( 26 );
            _clustalo_paramenters_tf.setText( "clustalo options" );
            _select_input_seqs_btn.addActionListener( this );
            _pnl.add( inputfile_pnl_1 );
            _pnl.add( inputfile_pnl_2 );
            _pnl.add( inputfile_pnl_3 );
            _pnl.add( inputfile_pnl_4 );
        }
        else {
            setTitle( "Phylogenetic Inference (from already aligned sequences) " );
            // Inputfile (MSA):
            final JPanel inputfile_pnl_1 = new JPanel();
            final JPanel inputfile_pnl_2 = new JPanel();
            inputfile_pnl_1.setLayout( new FlowLayout() );
            inputfile_pnl_2.setLayout( new FlowLayout() );
            inputfile_pnl_1.add( new JLabel( "Input MSA File:" ) );
            inputfile_pnl_1.add( _input_msa_file_tf = new JTextField() );
            inputfile_pnl_1.add( _select_input_msa_btn = new JButton( "Select Input File" ) );
            inputfile_pnl_2.add( new JLabel( "MSA: " ) );
            inputfile_pnl_2.add( new JLabel( "Number of Sequences:" ) );
            inputfile_pnl_2.add( _msa_size_tf = new JTextField() );
            inputfile_pnl_2.add( new JLabel( "Length:" ) );
            inputfile_pnl_2.add( _msa_length_tf = new JTextField() );
            inputfile_pnl_2.add( new JLabel( "Type:" ) );
            inputfile_pnl_2.add( _msa_type_tf = new JTextField() );
            _msa_length_tf.setColumns( 4 );
            _msa_size_tf.setColumns( 4 );
            _msa_type_tf.setColumns( 2 );
            _input_msa_file_tf.setColumns( 20 );
            _input_msa_file_tf.setEditable( false );
            _msa_length_tf.setEditable( false );
            _msa_size_tf.setEditable( false );
            _msa_type_tf.setEditable( false );
            _select_input_msa_btn.addActionListener( this );
            _pnl.add( inputfile_pnl_1 );
            _pnl.add( inputfile_pnl_2 );
        }
        //
        final JPanel inputfile_pnl_4 = new JPanel();
        inputfile_pnl_4.setLayout( new FlowLayout() );
        inputfile_pnl_4.add( new JLabel( "MSA Processing: " ) );
        inputfile_pnl_4.add( _execute_msa_processing_cb = new JCheckBox( "Process MSA" ) );
        inputfile_pnl_4.add( _msa_processing_remove_all_gap_columns_cb = new JCheckBox( "Remove all gap columns" ) );
        inputfile_pnl_4.add( new JLabel( "Max allowed gap ratio: " ) );
        inputfile_pnl_4.add( _msa_processing_max_allowed_gap_ratio_tf = new JTextField() );
        inputfile_pnl_4.add( new JLabel( "Min allowed non-gap sequence length: " ) );
        inputfile_pnl_4.add( _msa_processing_min_allowed_length_tf = new JTextField() );
        _msa_processing_max_allowed_gap_ratio_tf.setColumns( 4 );
        _msa_processing_min_allowed_length_tf.setColumns( 4 );
        final Border b = new LineBorder( Color.DARK_GRAY );
        inputfile_pnl_4.setBorder( b );
        _pnl.add( inputfile_pnl_4 );
        //
        // Distance calculation:
        // TODO if type==AA...
        final JPanel distance_calc_pnl_1 = new JPanel();
        distance_calc_pnl_1.setLayout( new FlowLayout() );
        distance_calc_pnl_1.add( new JLabel( "Distance calculation:" ) );
        distance_calc_pnl_1.add( _distance_calc_kimura_rb = new JRadioButton( "Kimura correction" ) );
        distance_calc_pnl_1.add( _distance_calc_poisson_rb = new JRadioButton( "Poisson" ) );
        distance_calc_pnl_1
                .add( _distance_calc_fract_dissimilarity_rb = new JRadioButton( "Fractional dissimilarity" ) );
        final ButtonGroup distance_calc_group_1 = new ButtonGroup();
        distance_calc_group_1.add( _distance_calc_kimura_rb );
        distance_calc_group_1.add( _distance_calc_poisson_rb );
        distance_calc_group_1.add( _distance_calc_fract_dissimilarity_rb );
        _pnl.add( distance_calc_pnl_1 );
        // Bootstrap resampling:
        final JPanel bootstrap_pnl = new JPanel();
        bootstrap_pnl.setLayout( new FlowLayout() );
        bootstrap_pnl.add( _bootstrap_cb = new JCheckBox( "Perform Bootstrap Resampling" ) );
        bootstrap_pnl.add( new JLabel( "Number of Bootstrap Samples:" ) );
        bootstrap_pnl.add( _bootstrap_tf = new JFormattedTextField( AptxUtil.createMaskFormatter( "###" ) ) );
        _bootstrap_tf.setColumns( 4 );
        // TODO see
        // http://download.oracle.com/javase/tutorial/uiswing/components/formattedtextfield.html
        // _bootstrap_tf.setColumns( 4 );
        bootstrap_pnl.add( new JLabel( "Random Seed:" ) );
        bootstrap_pnl.add( _random_seed_tf = new JTextField() );
        _random_seed_tf.setColumns( 4 );
        _pnl.add( bootstrap_pnl );
        final JPanel launch_pnl = new JPanel();
        launch_pnl.setLayout( new FlowLayout() );
        _launch_btn = new JButton( "Go!" );
        _launch_btn.addActionListener( this );
        launch_pnl.add( _launch_btn );
        _cancel_btn = new JButton( "Cancel" );
        _cancel_btn.addActionListener( this );
        launch_pnl.add( _cancel_btn );
        _pnl.add( launch_pnl );
        initializeValues( from_unaligned_seqs );
        pack();
        setLocationRelativeTo( getParentFrame() );
        setResizable( false );
    }

    @Override
    public void actionPerformed( final ActionEvent e ) {
        if ( e.getSource() == _select_input_msa_btn ) {
            readInputFile();
        }
        else if ( e.getSource() == _select_input_seqs_btn ) {
            readInputSeqsFile();
        }
        else if ( e.getSource() == _launch_btn ) {
            launch();
        }
        else if ( e.getSource() == _cancel_btn ) {
            cancel();
        }
    }

    public void activate() {
        setVisible( true );
    }

    private MainFrameApplication getParentFrame() {
        return _parent_frame;
    }

    public PhylogeneticInferenceOptions getPhylogeneticInferenceOptions() {
        return _opts;
    }

    public int getValue() {
        return _value;
    }

    private void initializeValues( final boolean from_unaligned_seqs ) {
        _value = JOptionPane.CANCEL_OPTION;
        if ( from_unaligned_seqs ) {
            updateSeqsItems();
        }
        else {
            updateMsaItems();
        }
        updateMsaProcessingItem();
        updateDistanceCalcMethod();
        _bootstrap_tf.setText( getPhylogeneticInferenceOptions().getBootstrapSamples() + "" );
        _random_seed_tf.setText( getPhylogeneticInferenceOptions().getRandomNumberGeneratorSeed() + "" );
    }

    private void launch() {
        processPerformBootstrapResampling();
        if ( _bootstrap_cb.isSelected() ) {
            processBootstrapSamplesNumber();
            processRandomNumberGeneratorSeed();
        }
        if ( true ) {
            //TODO
            processMsaProcessing();
        }
        processDistanceCalcMethod();
        processMsaPrgParameters();
        setVisible( false );
        _value = JOptionPane.OK_OPTION;
    }

    private void cancel() {
        setVisible( false );
        _value = JOptionPane.CANCEL_OPTION;
    }

    private void processBootstrapSamplesNumber() {
        int bootstrap_samples = 0;
        try {
            bootstrap_samples = Integer.parseInt( _bootstrap_tf.getText().trim() );
        }
        catch ( final NumberFormatException e ) {
            // JOptionPane.showMessageDialog( this, "Could not parse number of bootstrap resamplings from: " +  _bootstrap_tf.getText().trim(), "User Error", JOptionPane.ERROR_MESSAGE );
            return;
        }
        if ( bootstrap_samples >= 0 ) {
            getPhylogeneticInferenceOptions().setBootstrapSamples( bootstrap_samples );
        }
    }

    private void processRandomNumberGeneratorSeed() {
        long seed = PhylogeneticInferenceOptions.RANDOM_NUMBER_SEED_DEFAULT;
        try {
            seed = Long.parseLong( _random_seed_tf.getText().trim() );
        }
        catch ( final NumberFormatException e ) {
            return;
        }
        getPhylogeneticInferenceOptions().setRandomNumberGeneratorSeed( seed );
    }

    private void processMsaProcessing() {
        getPhylogeneticInferenceOptions().setExecuteMsaProcessing( _execute_msa_processing_cb.isSelected() );
        getPhylogeneticInferenceOptions()
                .setMsaProcessingRemoveAllGapColumns( _msa_processing_remove_all_gap_columns_cb.isSelected() );
        int min_length = -1;
        try {
            min_length = Integer.parseInt( _msa_processing_min_allowed_length_tf.getText().trim() );
        }
        catch ( final NumberFormatException e ) {
            min_length = -1;
        }
        if ( min_length > 0 ) {
            getPhylogeneticInferenceOptions().setMsaProcessingMinAllowedLength( min_length );
        }
        double msa_processing_max_allowed_gap_ratio = -1.0;
        try {
            msa_processing_max_allowed_gap_ratio = Double.parseDouble( _msa_processing_max_allowed_gap_ratio_tf
                    .getText().trim() );
        }
        catch ( final NumberFormatException e ) {
            msa_processing_max_allowed_gap_ratio = -1.0;
        }
        if ( ( msa_processing_max_allowed_gap_ratio >= 0.0 ) && ( msa_processing_max_allowed_gap_ratio <= 1.0 ) ) {
            getPhylogeneticInferenceOptions().setMsaProcessingMaxAllowedGapRatio( msa_processing_max_allowed_gap_ratio );
        }
    }

    private void processDistanceCalcMethod() {
        if ( ( _distance_calc_kimura_rb != null ) && _distance_calc_kimura_rb.isSelected() ) {
            getPhylogeneticInferenceOptions().setPwdDistanceMethod( PWD_DISTANCE_METHOD.KIMURA_DISTANCE );
        }
        else if ( ( _distance_calc_poisson_rb != null ) && _distance_calc_poisson_rb.isSelected() ) {
            getPhylogeneticInferenceOptions().setPwdDistanceMethod( PWD_DISTANCE_METHOD.POISSON_DISTANCE );
        }
        else if ( ( _distance_calc_fract_dissimilarity_rb != null )
                && _distance_calc_fract_dissimilarity_rb.isSelected() ) {
            getPhylogeneticInferenceOptions().setPwdDistanceMethod( PWD_DISTANCE_METHOD.FRACTIONAL_DISSIMILARITY );
        }
    }

    private void processPerformBootstrapResampling() {
        getPhylogeneticInferenceOptions().setPerformBootstrapResampling( _bootstrap_cb.isSelected() );
    }

    private void processMsaPrgParameters() {
        if ( _mafft_paramenters_tf != null ) {
            getPhylogeneticInferenceOptions().setMsaPrgParameters( _mafft_paramenters_tf.getText() );
        }
    }

    private void readInputFile() {
        getParentFrame().readMsaFromFile();
        updateMsaItems();
    }

    private void readInputSeqsFile() {
        getParentFrame().readSeqsFromFile();
        updateSeqsItems();
    }

    private void updateDistanceCalcMethod() {
        switch ( getPhylogeneticInferenceOptions().getPwdDistanceMethod() ) {
            case KIMURA_DISTANCE:
                _distance_calc_kimura_rb.setSelected( true );
                break;
            case POISSON_DISTANCE:
                _distance_calc_poisson_rb.setSelected( true );
                break;
            case FRACTIONAL_DISSIMILARITY:
                _distance_calc_fract_dissimilarity_rb.setSelected( true );
                break;
            default:
                throw new RuntimeException( "invalid distance calc method" );
        }
    }

    private void updateMsaProcessingItem() {
        _execute_msa_processing_cb.setSelected( getPhylogeneticInferenceOptions().isExecuteMsaProcessing() );
        _msa_processing_remove_all_gap_columns_cb.setSelected( getPhylogeneticInferenceOptions()
                .isMsaProcessingRemoveAllGapColumns() );
        if ( _opts.getMsaProcessingMaxAllowedGapRatio() > 0 ) {
            _msa_processing_max_allowed_gap_ratio_tf.setText( _opts.getMsaProcessingMaxAllowedGapRatio() + "" );
        }
        if ( _opts.getMsaProcessingMinAllowedLength() > 0 ) {
            _msa_processing_min_allowed_length_tf.setText( _opts.getMsaProcessingMinAllowedLength() + "" );
        }
    }

    private void updateMsaItems() {
        if ( getParentFrame().getMsa() != null ) {
            _input_msa_file_tf.setText( getParentFrame().getMsaFile().toString() );
            _msa_length_tf.setText( getParentFrame().getMsa().getLength() + "" );
            _msa_size_tf.setText( getParentFrame().getMsa().getNumberOfSequences() + "" );
            _msa_type_tf.setText( getParentFrame().getMsa().getType() + "" );
            _input_msa_file_tf.setEnabled( true );
            _msa_length_tf.setEnabled( true );
            _msa_size_tf.setEnabled( true );
            _msa_type_tf.setEnabled( true );
            _launch_btn.setEnabled( true );
        }
        else {
            _input_msa_file_tf.setText( "" );
            _msa_length_tf.setText( "" );
            _msa_size_tf.setText( "" );
            _msa_type_tf.setText( "" );
            _input_msa_file_tf.setEnabled( false );
            _msa_length_tf.setEnabled( false );
            _msa_size_tf.setEnabled( false );
            _msa_type_tf.setEnabled( false );
            _launch_btn.setEnabled( false );
        }
    }

    private void updateSeqsItems() {
        if ( getParentFrame().getSeqs() != null ) {
            final DescriptiveStatistics stats = calcSequenceStats( getParentFrame().getSeqs() );
            _input_seqs_tf.setText( getParentFrame().getSeqsFile().toString() );
            _input_seqs_median_length_tf.setText( ( int ) stats.median() + "" );
            _input_seqs_min_length_tf.setText( ( int ) stats.getMin() + "" );
            _input_seqs_max_length_tf.setText( ( int ) stats.getMax() + "" );
            _input_seqs_number_tf.setText( getParentFrame().getSeqs().size() + "" );
            _input_seqs_type_tf.setText( getParentFrame().getSeqs().get( 0 ).getType() + "" );
            _input_seqs_tf.setEnabled( true );
            _input_seqs_median_length_tf.setEnabled( true );
            _input_seqs_min_length_tf.setEnabled( true );
            _input_seqs_max_length_tf.setEnabled( true );
            _input_seqs_number_tf.setEnabled( true );
            _input_seqs_type_tf.setEnabled( true );
            _launch_btn.setEnabled( true );
        }
        else {
            _input_seqs_tf.setText( "" );
            _input_seqs_median_length_tf.setText( "" );
            _input_seqs_min_length_tf.setText( "" );
            _input_seqs_max_length_tf.setText( "" );
            _input_seqs_number_tf.setText( "" );
            _input_seqs_type_tf.setText( "" );
            _input_seqs_tf.setEnabled( false );
            _input_seqs_median_length_tf.setEnabled( false );
            _input_seqs_min_length_tf.setEnabled( false );
            _input_seqs_max_length_tf.setEnabled( false );
            _input_seqs_number_tf.setEnabled( false );
            _input_seqs_type_tf.setEnabled( false );
            _launch_btn.setEnabled( false );
        }
    }

    DescriptiveStatistics calcSequenceStats( final List<Sequence> seqs ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final Sequence s : seqs ) {
            stats.addValue( s.getLength() );
        }
        return stats;
    }
}
