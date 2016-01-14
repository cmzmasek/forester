
package org.forester.io.parsers;

import java.io.IOException;

import org.forester.phylogeny.Phylogeny;

public interface IteratingPhylogenyParser {

    public void reset() throws IOException;

    public Phylogeny next() throws IOException;

    public boolean hasNext();

    public void setSource( final Object o ) throws IOException;
}
