package org.broadinstitute.sting.gatk.arguments;

import org.broadinstitute.sting.commandline.Argument;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/13/12
 * Time: 10:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class IonArgumentCollection extends GATKArgumentCollection {

    @Argument(fullName = "bypassFlowAlign", shortName = "nofa", doc = "Bypasses Flow Alignment calculation.", required = false)
    public boolean bypassFlowAlign = false;

}
