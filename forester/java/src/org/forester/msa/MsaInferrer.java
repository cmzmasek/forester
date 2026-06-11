// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.msa;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.forester.sequence.MolecularSequence;
import org.forester.util.SystemCommandExecutor;

public abstract class MsaInferrer {

    public abstract String getErrorDescription();

    public abstract int getExitCode();

    public static boolean isInstalled( final String path_to_prg ) {
        return SystemCommandExecutor.isExecuteableFile( new File( path_to_prg ) );
    }

    @Override
    public Object clone() {
        throw new NoSuchMethodError();
    }

    public abstract Msa infer( File path_to_input_seqs, List<String> opts ) throws IOException, InterruptedException;

    public abstract Msa infer( final List<MolecularSequence> seqs, final List<String> opts ) throws IOException,
    InterruptedException;
}
