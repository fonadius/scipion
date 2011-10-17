/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JDialogXmippFilesList.java
 *
 * Created on Aug 2, 2011, 5:56:47 PM
 */
package browser.filebrowsers;

import browser.COMMAND_PARAMETERS;
import browser.LABELS;
import browser.imageitems.listitems.FileItem;
import browser.windows.ImagesWindowFactory;
import ij.IJ;
import java.awt.BorderLayout;
import java.io.OutputStreamWriter;
import java.net.InetAddress;
import java.net.Socket;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

/**
 *
 * @author Juanjo Vega
 */
public class JDialogXmippFilesList extends javax.swing.JFrame {

    final static String SOT = "__STARTED__";
    final static String EOT = "__END__";
    String seltype;
    protected JPanelXmippBrowser panelXmippBrowser;
    int port;

    public JDialogXmippFilesList(String directory, int port) {
        this(directory, port, false, COMMAND_PARAMETERS.SELECTION_TYPE_ANY, "");
    }

    public JDialogXmippFilesList(final String directory, int port,
            final boolean singleSelection, String seltype, final String expression) {
        super();

        this.port = port;
        this.seltype = seltype;

        // Tells listener that it has been started.
        send(new Object[]{SOT}, false);

        if (seltype.compareTo(COMMAND_PARAMETERS.SELECTION_TYPE_FILE) == 0) {
            setTitle(LABELS.TITLE_XMIPP_FILE_SELECTOR_FILE);
        } else if (seltype.compareTo(COMMAND_PARAMETERS.SELECTION_TYPE_DIR) == 0) {
            setTitle(LABELS.TITLE_XMIPP_FILE_SELECTOR_DIR);
        } else {
            setTitle(LABELS.TITLE_XMIPP_FILE_SELECTOR_ANY);
        }

        initComponents();

        ImagesWindowFactory.setConvenientSize(this);
        this.setLocationRelativeTo(null);

        setPanel(directory, expression, singleSelection);
    }

    protected void setPanel(final String directory, final String expression, final boolean singleSelection) {
        ImagesWindowFactory.blockGUI(getRootPane(), "Building list...");

        Thread t = new Thread(new Runnable() {

            public void run() {
                panelXmippBrowser = new JPanelXmippBrowser(directory, expression);
                panelXmippBrowser.setSingleSelection(singleSelection);

                add(panelXmippBrowser, BorderLayout.CENTER);

                pack();

                ImagesWindowFactory.releaseGUI(getRootPane());
            }
        });

        t.start();

        try {
            long time = System.currentTimeMillis();

            t.join();

            time = System.currentTimeMillis() - time;
            SimpleDateFormat dateFormatter = new SimpleDateFormat("'List building time:' mm:ss:S");
            System.out.println(dateFormatter.format(new Date(time)));
        } catch (Exception ex) {
        }
    }

    protected void button1Clicked() {
        if (sendSelectedFiles()) {
            dispose();
        }

//        DEBUG.printMessage("Exiting...");
        System.exit(0);
    }

    protected void button2Clicked() {
        cancel();
    }

    protected void cancel() {
        send(null, true);
        dispose();

//        DEBUG.printMessage("Exiting...");
        System.exit(0);
    }

    protected boolean sendSelectedFiles() {
        Object objs[] = panelXmippBrowser.getSelectedValues();
        ArrayList<FileItem> list = new ArrayList<FileItem>();

        for (int i = 0; i < objs.length; i++) {
            FileItem item = (FileItem) objs[i];

            if (acceptFile(item)) {
                list.add(item);
            }
        }

        return send(list.toArray(), true);
    }

    boolean acceptFile(FileItem item) {
        if (seltype.compareTo(COMMAND_PARAMETERS.SELECTION_TYPE_ANY) == 0) {
            return true;
        } else if (item.isDirectory()) {
            return seltype.compareTo(COMMAND_PARAMETERS.SELECTION_TYPE_DIR) == 0;
        } else {
            return seltype.compareTo(COMMAND_PARAMETERS.SELECTION_TYPE_FILE) == 0;
        }
    }

    final protected boolean send(Object items[], boolean end) {
        try {
            Socket socket = new Socket(InetAddress.getByName("127.0.0.1"), port);

            // Get streams.
            OutputStreamWriter output = new OutputStreamWriter(socket.getOutputStream());
            output.flush();

            if (items != null) {
                for (int i = 0; i < items.length; i++) {
                    output.write(items[i].toString() + "\n");
                }
            }

            if (end) {
                output.write(EOT + "\n");
            }

            output.flush();

            // Closes connection.
            output.close();
            socket.close();

            return true;
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }

        return false;
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jpButtons = new javax.swing.JPanel();
        jbOk = new javax.swing.JButton();
        jbCancel = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
        });

        jbOk.setText(LABELS.BUTTON_OK);
        jbOk.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbOkActionPerformed(evt);
            }
        });
        jpButtons.add(jbOk);

        jbCancel.setText(LABELS.BUTTON_CANCEL);
        jbCancel.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbCancelActionPerformed(evt);
            }
        });
        jpButtons.add(jbCancel);

        getContentPane().add(jpButtons, java.awt.BorderLayout.SOUTH);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jbOkActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbOkActionPerformed
        button1Clicked();
}//GEN-LAST:event_jbOkActionPerformed

    private void jbCancelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbCancelActionPerformed
        button2Clicked();
}//GEN-LAST:event_jbCancelActionPerformed

    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        cancel();
    }//GEN-LAST:event_formWindowClosing
    // Variables declaration - do not modify//GEN-BEGIN:variables
    protected javax.swing.JButton jbCancel;
    protected javax.swing.JButton jbOk;
    protected javax.swing.JPanel jpButtons;
    // End of variables declaration//GEN-END:variables
}
