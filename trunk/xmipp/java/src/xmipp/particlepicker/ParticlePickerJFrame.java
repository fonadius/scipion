package xmipp.particlepicker;

import ij.CommandListener;
import ij.Executer;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.frame.Recorder;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import xmipp.ij.commons.Tool;
import xmipp.ij.commons.XmippIJUtil;
import xmipp.particlepicker.tiltpair.gui.TiltPairParticlesJDialog;
import xmipp.particlepicker.training.gui.TrainingPickerJFrame;
import xmipp.particlepicker.training.model.FamilyState;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.utils.ColorIcon;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;

public abstract class ParticlePickerJFrame extends JFrame implements ActionListener
{

	protected ParticlesJDialog particlesdialog;

	protected JMenuItem ijmi;
	protected JCheckBox circlechb;
	protected JCheckBox rectanglechb;
	protected JFormattedTextField sizetf;
	protected JCheckBox centerchb;
	protected JPanel symbolpn;
	protected JMenuItem savemi;
	protected JMenuItem hcontentsmi;
	protected JMenuItem pmi;
	protected JMenuItem importfpmi;
	protected JMenuItem exportmi;
	protected JMenu filtersmn;
	protected String activefilter;
	protected JSlider sizesl;
	protected JPanel sizepn;
	private String command;
	private List<JCheckBoxMenuItem> mifilters;
	protected JMenu filemn;
	protected JMenuItem importffmi;
	protected JButton colorbt;
	protected Color color;
	protected JPanel colorpn;
	protected JButton resetbt;

	private JMenuItem exitmi;

	public ParticlePickerJFrame(ParticlePicker picker)
	{

		addWindowListener(new WindowAdapter()
		{
			public void windowClosing(WindowEvent winEvt)
			{
				if (getParticlePicker().isChanged() && XmippDialog.showQuestion(ParticlePickerJFrame.this, "Save changes before closing?"))
					saveChanges();
				System.exit(0);
			}
		});

		Recorder.record = true;

		// detecting if a command is thrown by ImageJ
		Executer.addCommandListener(new CommandListener()
		{
			public String commandExecuting(String command)
			{
				ParticlePickerJFrame.this.command = command;
				return command;

			}
		});
		ImagePlus.addImageListener(new ImageListener()
		{

			@Override
			public void imageUpdated(ImagePlus arg0)
			{
				if (command != null)
				{
					String options = "";
					if (Recorder.getCommandOptions() != null)
						options = Recorder.getCommandOptions();
					if (!getParticlePicker().isFilterSelected(command))
						getParticlePicker().addFilter(command, options);
					command = null;
				}
			}

			@Override
			public void imageOpened(ImagePlus arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void imageClosed(ImagePlus arg0)
			{
				// TODO Auto-generated method stub

			}
		});
		filemn = new JMenu("File");
		savemi = new JMenuItem("Save", XmippResource.getIcon("save.gif"));
		savemi.setMnemonic('S');
		savemi.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_DOWN_MASK));
		savemi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				saveChanges();
				showMessage("Data saved successfully");
				((JMenuItem) e.getSource()).setEnabled(false);
			}
		});
		filemn.add(savemi);
		importfpmi = new JMenuItem("Import from Project...", XmippResource.getIcon("import_wiz.gif"));
		filemn.add(importfpmi);
		importfpmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				new ImportParticlesFromProjectJDialog(ParticlePickerJFrame.this, true);
			}
		});

		importffmi = new JMenuItem("Import from Files...", XmippResource.getIcon("import.gif"));
		filemn.add(importffmi);
		importffmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				displayImportDialog();

			}
		});

		exportmi = new JMenuItem("Export Particles...", XmippResource.getIcon("export_wiz.gif"));
		filemn.add(exportmi);
		exportmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				XmippFileChooser fc = new XmippFileChooser();
				int returnVal = fc.showOpenDialog(ParticlePickerJFrame.this);

				try
				{
					if (returnVal == XmippFileChooser.APPROVE_OPTION)
					{
						File file = fc.getSelectedFile();
						getParticlePicker().exportParticles(getFamily(), file.getAbsolutePath());
						showMessage("Export successful");
					}
				}
				catch (Exception ex)
				{
					showException(ex);
				}
			}
		});
		exitmi = new JMenuItem("Exit");
		exitmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent arg0)
			{

			}
		});
		filemn.add(exitmi);

		ijmi = new JMenuItem("ImageJ", XmippResource.getIcon("ij.gif"));
		ijmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				XmippIJUtil.showImageJ(Tool.PICKER);
			}
		});

		hcontentsmi = new JMenuItem("Online help", XmippResource.getIcon("online_help.gif"));
		hcontentsmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				try
				{
					XmippWindowUtil.openURI("http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParticlePicker");
				}
				catch (Exception ex)
				{
					showException(ex);
					// JOptionPane.showMessageDialog(ParticlePickerJFrame.this,
					// ex.getMessage());
				}
			}
		});
		pmi = new JMenuItem("Particles", XmippResource.getIcon("table_view.gif"));
		pmi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				loadParticles();
			}
		});

		mifilters = new ArrayList<JCheckBoxMenuItem>();
		filtersmn = new JMenu("Filters");
		filtersmn.addMenuListener(new MenuListener()
		{

			@Override
			public void menuCanceled(MenuEvent arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void menuDeselected(MenuEvent arg0)
			{
				// TODO Auto-generated method stub

			}

			@Override
			public void menuSelected(MenuEvent arg0)
			{
				for (JCheckBoxMenuItem mi : mifilters)
					mi.setSelected(getParticlePicker().isFilterSelected(mi.getText()));

			}
		});

		addFilterMenuItem("Smooth Filter", true, picker);
		addFilterMenuItem("Bandpass Filter...", true, picker);

		JCheckBoxMenuItem admi = addFilterMenuItem("Anisotropic Diffusion...", false, picker);
		admi.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				activefilter = "8-bit";
				IJ.run(activefilter);
				activefilter = ((JCheckBoxMenuItem) e.getSource()).getText();
				IJ.run(activefilter);
			}
		});
		addFilterMenuItem("Mean Shift", true, picker);
		addFilterMenuItem("Subtract Background...", true, picker);
		addFilterMenuItem("Gaussian Blur...", true, picker);
		addFilterMenuItem("Brightness/Contrast...", true, picker);

		resetbt = XmippWindowUtil.getTextButton("Reset", new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				resetMicrograph();
			}
		});

	}

	//
	// protected void smooth()
	// {
	// try
	// {
	// ImageGeneric Iaux =
	// XmippImageConverter.convertToImageGeneric(IJ.getImage());
	// Iaux.convert2Datatype(ImageGeneric.UChar);
	// ImageGeneric Ismooth = new ImageGeneric(ImageGeneric.UChar);
	// Ismooth.resize(Iaux.getXDim(), Iaux.getYDim());
	// Iaux.smooth(Ismooth);
	// ImagePlus imp = XmippImageConverter.convertToImagePlus(Ismooth);
	// }
	// catch (Exception e1)
	// {
	// // TODO Auto-generated catch block
	// e1.printStackTrace();
	// }
	// }

	protected abstract void resetMicrograph();

	protected void enableEdition(boolean enable)
	{
		importffmi.setEnabled(enable);
		importfpmi.setEnabled(enable);
		savemi.setEnabled(enable);
		sizesl.setEnabled(enable);
		colorbt.setEnabled(enable);
		resetbt.setEnabled(enable);
	}

	protected abstract void displayImportDialog();

	private JCheckBoxMenuItem addFilterMenuItem(String command, boolean defaultlistener, ParticlePicker picker)
	{
		JCheckBoxMenuItem mi = new JCheckBoxMenuItem(command);
		mifilters.add(mi);
		mi.setSelected(picker.isFilterSelected(command));
		if (defaultlistener)
			mi.addActionListener(this);
		filtersmn.add(mi);
		return mi;
	}

	@Override
	public void actionPerformed(ActionEvent e)
	{
		try
		{
			JCheckBoxMenuItem item = (JCheckBoxMenuItem) e.getSource();
			activefilter = item.getText();
			if (item.isSelected())// filter added, will be registered by picker
									// with options if needed
				if (activefilter.equals("Smooth Filter"))
				{
					getParticlePicker().addFilter("Smooth Filter", "xmipp");
					reloadImage();
				}
				else
				{
					for(int i = 0; i < WindowManager.getImageCount(); i ++)
						IJ.run(WindowManager.getImage(i), activefilter, "");
				}
			else
			{
				// filter removed
				getParticlePicker().removeFilter(activefilter);
				reloadImage();
			}

			setChanged(true);
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
			showException(ex);
		}

	}
	
	protected abstract void reloadImage();


	protected abstract void saveChanges();

	public int getSide(int size)
	{
		return 100;
	}

	public abstract Family getFamily();

	public abstract ParticlePickerCanvas getCanvas();

	public void loadParticles()
	{
		try
		{
			if (particlesdialog == null)
				if (ParticlePickerJFrame.this instanceof TrainingPickerJFrame)
					particlesdialog = new ParticlesJDialog(ParticlePickerJFrame.this);
				else
					particlesdialog = new TiltPairParticlesJDialog(ParticlePickerJFrame.this);
			else
			{

				particlesdialog.loadParticles(true);
				particlesdialog.setVisible(true);
			}
		}
		catch (Exception ex)
		{
			showException(ex);
			if (particlesdialog != null)
				particlesdialog.close();
			particlesdialog = null;
		}
	}

	public void updateMicrographsModel()
	{
		if (particlesdialog != null)
			loadParticles();
	}

	public ParticlesJDialog getParticlesJDialog()
	{
		return particlesdialog;
	}

	public abstract Micrograph getMicrograph();

	public abstract List<? extends TrainingParticle> getParticles();

	public boolean isPickingAvailable(MouseEvent e)
	{
		if (getCanvas().getTool() != Tool.PICKER)
			return false;
		if (SwingUtilities.isRightMouseButton(e))
			return false;
		if (getParticlePicker().getMode() == FamilyState.ReadOnly)
			return false;
		return true;
	}

	protected void initSymbolPane()
	{

		symbolpn = new JPanel(new FlowLayout(FlowLayout.LEFT));
		symbolpn.add(new JLabel("Symbol:"));
		// symbolpn.setBorder(BorderFactory.createTitledBorder("Symbol"));
		ShapeItemListener shapelistener = new ShapeItemListener();

		circlechb = new JCheckBox(Shape.Circle.toString());
		circlechb.setSelected(true);
		circlechb.addItemListener(shapelistener);

		rectanglechb = new JCheckBox(Shape.Rectangle.toString());
		rectanglechb.setSelected(true);
		rectanglechb.addItemListener(shapelistener);

		centerchb = new JCheckBox(Shape.Center.toString());
		centerchb.setSelected(true);
		centerchb.addItemListener(shapelistener);

		symbolpn.add(circlechb);
		symbolpn.add(rectanglechb);
		symbolpn.add(centerchb);
	}

	class ShapeItemListener implements ItemListener
	{
		@Override
		public void itemStateChanged(ItemEvent e)
		{
			changeShapes();
		}
	}

	public abstract void changeShapes();

	public boolean isShapeSelected(Shape shape)
	{
		switch (shape)
		{
		case Rectangle:
			return rectanglechb.isSelected();
		case Circle:
			return circlechb.isSelected();
		case Center:
			return centerchb.isSelected();
			// case OnlyLast:
			// return onlylastchb.isSelected();
		}
		return false;
	}

	public abstract ParticlePicker getParticlePicker();

	public abstract void setChanged(boolean changed);

	protected void initColorPane()
	{
		colorpn = new JPanel();
		color = getFamily().getColor();
		colorpn.add(new JLabel("Color:"));
		colorbt = new JButton();
		colorbt.setIcon(new ColorIcon(color));
		colorbt.setBorderPainted(false);
		colorpn.add(colorbt);
	}

	protected void initSizePane()
	{
		sizepn = new JPanel();

		int size = getFamily().getSize();
		sizepn.add(new JLabel("Size:"));
		sizesl = new JSlider(0, 1000, size);
		sizesl.setPaintTicks(true);
		sizesl.setMajorTickSpacing(100);
		int height = (int) sizesl.getPreferredSize().getHeight();

		sizesl.setPreferredSize(new Dimension(100, height));
		sizepn.add(sizesl);

		sizetf = new JFormattedTextField(NumberFormat.getIntegerInstance());
		sizetf.setColumns(3);
		sizetf.setText(Integer.toString(size));
		sizepn.add(sizetf);
		sizetf.addActionListener(new ActionListener()
		{

			@Override
			public void actionPerformed(ActionEvent e)
			{
				int size = ((Number) sizetf.getValue()).intValue();
				switchSize(size);
			}
		});

		sizesl.addChangeListener(new ChangeListener()
		{

			@Override
			public void stateChanged(ChangeEvent e)
			{
				int size = sizesl.getValue();
				switchSize(size);
			}
		});

	}

	public void switchSize(int size)
	{
		sizetf.setText(Integer.toString(size));
		sizesl.setValue(size);
		getCanvas().repaint();
		getFamily().setSize(size);
		if (particlesdialog != null)
		{
			for (TrainingParticle p : getParticles())
				p.resetParticleCanvas();
			loadParticles();
		}
		setChanged(true);
	}

	protected void importParticlesXmipp30(String path)
	{
		getParticlePicker().importParticlesXmipp30(getFamily(), path);
		setChanged(true);
		getCanvas().repaint();
		updateMicrographsModel();
	}

	public void importParticlesXmipp24(String projectdir)
	{
		getParticlePicker().importParticlesFromXmipp24Project(getFamily(), projectdir);
		setChanged(true);
		getCanvas().repaint();
		updateMicrographsModel();

	}

	public void importParticlesEman(String path)
	{
		throw new UnsupportedOperationException(XmippMessage.getNotImplementedYetMsg());

	}

	/** Shortcut function to show messages */
	private boolean showMessage(String message)
	{
		return XmippDialog.showInfo(this, message);
	}

	private boolean showException(Exception e)
	{
		return XmippDialog.showException(this, e);
	}

}
