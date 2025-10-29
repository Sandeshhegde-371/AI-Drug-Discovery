import { HelpCircle, FileText, Cpu, BarChart, Info } from 'lucide-react';
import { Card } from './ui/card';
import {
  Accordion,
  AccordionContent,
  AccordionItem,
  AccordionTrigger,
} from './ui/accordion';

export function InstructionsPanel() {
  return (
    <Card className="p-6 bg-white/70 backdrop-blur-sm border-purple-100 shadow-lg shadow-purple-100/50">
      <div className="flex items-center gap-2 mb-4">
        <div className="w-10 h-10 bg-gradient-to-br from-amber-500 to-orange-500 rounded-lg flex items-center justify-center">
          <HelpCircle className="w-5 h-5 text-white" />
        </div>
        <div>
          <h3 className="text-gray-900">Quick Guide</h3>
          <p className="text-sm text-gray-600">Everything you need to know</p>
        </div>
      </div>

      <Accordion type="single" collapsible className="space-y-2">
        <AccordionItem value="getting-started" className="border border-purple-100 rounded-lg px-4 bg-white/50">
          <AccordionTrigger className="text-left hover:no-underline py-3">
            <div className="flex items-center gap-2">
              <Info className="w-4 h-4 text-purple-600" />
              <span>Getting Started</span>
            </div>
          </AccordionTrigger>
          <AccordionContent className="pb-4">
            <div className="space-y-3 text-sm text-gray-700">
              <p>Follow these simple steps to analyze your data:</p>
              <ol className="space-y-2 ml-4 list-decimal">
                <li>Upload your CSV file with any column structure</li>
                <li>Select the AI pipeline that matches your goal</li>
                <li>Configure parameters to fine-tune the analysis</li>
                <li>Run the pipeline and download your results</li>
              </ol>
              <div className="mt-3 p-3 bg-purple-50 border border-purple-200 rounded-lg">
                <p className="text-purple-900 text-xs">
                  💡 <strong>Tip:</strong> Start with the sample CSV to familiarize yourself with the platform
                </p>
              </div>
            </div>
          </AccordionContent>
        </AccordionItem>

        <AccordionItem value="file-format" className="border border-purple-100 rounded-lg px-4 bg-white/50">
          <AccordionTrigger className="text-left hover:no-underline py-3">
            <div className="flex items-center gap-2">
              <FileText className="w-4 h-4 text-purple-600" />
              <span>CSV File Format</span>
            </div>
          </AccordionTrigger>
          <AccordionContent className="pb-4">
            <div className="space-y-3 text-sm text-gray-700">
              <p>Your CSV file can contain any columns with any data structure. Example format:</p>
              <div className="bg-gray-900 text-green-400 p-3 rounded-lg border border-gray-700 font-mono text-xs overflow-x-auto">
                <div>compound_id,structure,weight,property_1,...</div>
                <div className="text-gray-500">COMP001,CC(=O)OC1=CC=CC=C1,180.16,1.19,...</div>
              </div>
              <ul className="space-y-2 ml-4">
                <li className="flex items-start gap-2">
                  <span className="text-purple-600 mt-1">•</span>
                  <span>Any column names are accepted (molecular data, properties, numeric values, etc.)</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-purple-600 mt-1">•</span>
                  <span>At least one header row and one data row required</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-purple-600 mt-1">•</span>
                  <span>Use comma separation only (no tabs or semicolons)</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-purple-600 mt-1">•</span>
                  <span>All columns preserved in the output</span>
                </li>
              </ul>
            </div>
          </AccordionContent>
        </AccordionItem>

        <AccordionItem value="pipelines" className="border border-purple-100 rounded-lg px-4 bg-white/50">
          <AccordionTrigger className="text-left hover:no-underline py-3">
            <div className="flex items-center gap-2">
              <Cpu className="w-4 h-4 text-purple-600" />
              <span>Understanding Pipelines</span>
            </div>
          </AccordionTrigger>
          <AccordionContent className="pb-4">
            <div className="space-y-4 text-sm text-gray-700">
              <div className="p-3 bg-purple-50 border-l-4 border-purple-600 rounded">
                <h5 className="text-purple-900 mb-1">Generative Screening</h5>
                <p className="text-gray-700">Creates new molecular structures based on your seed data using deep learning models. Ideal for expanding chemical libraries and discovering novel candidates.</p>
              </div>

              <div className="p-3 bg-blue-50 border-l-4 border-blue-600 rounded">
                <h5 className="text-blue-900 mb-1">Virtual Screening</h5>
                <p className="text-gray-700">Rapidly evaluates large compound libraries to identify promising candidates. Uses docking simulations and binding affinity predictions.</p>
              </div>

              <div className="p-3 bg-green-50 border-l-4 border-green-600 rounded">
                <h5 className="text-green-900 mb-1">Predict Properties</h5>
                <p className="text-gray-700">Estimates physicochemical and biological properties including activity, solubility, toxicity, and ADMET characteristics using ML models.</p>
              </div>

              <div className="p-3 bg-orange-50 border-l-4 border-orange-600 rounded">
                <h5 className="text-orange-900 mb-1">Prioritize Candidates</h5>
                <p className="text-gray-700">Ranks molecules using multi-criteria decision analysis to optimize which compounds to test in the lab first, saving time and resources.</p>
              </div>
            </div>
          </AccordionContent>
        </AccordionItem>

        <AccordionItem value="output" className="border border-purple-100 rounded-lg px-4 bg-white/50">
          <AccordionTrigger className="text-left hover:no-underline py-3">
            <div className="flex items-center gap-2">
              <BarChart className="w-4 h-4 text-purple-600" />
              <span>Understanding Results</span>
            </div>
          </AccordionTrigger>
          <AccordionContent className="pb-4">
            <div className="space-y-3 text-sm text-gray-700">
              <p>Each pipeline adds specific columns to your output:</p>
              <div className="space-y-2">
                <div className="p-2 bg-gray-50 rounded border border-gray-200">
                  <strong className="text-purple-600">Score:</strong> Confidence metric (0-1 scale)
                </div>
                <div className="p-2 bg-gray-50 rounded border border-gray-200">
                  <strong className="text-purple-600">High Confidence:</strong> Score {'>'} 0.8
                </div>
                <div className="p-2 bg-gray-50 rounded border border-gray-200">
                  <strong className="text-purple-600">Medium Confidence:</strong> Score 0.6-0.8
                </div>
                <div className="p-2 bg-gray-50 rounded border border-gray-200">
                  <strong className="text-purple-600">Low Confidence:</strong> Score {'<'} 0.6
                </div>
              </div>
              <div className="mt-3 p-3 bg-amber-50 border border-amber-200 rounded-lg">
                <p className="text-amber-900 text-xs">
                  ⚠️ <strong>Important:</strong> Always validate AI predictions with experimental data before making research decisions
                </p>
              </div>
            </div>
          </AccordionContent>
        </AccordionItem>
      </Accordion>

      <div className="mt-6 p-4 bg-gradient-to-br from-purple-50 to-indigo-50 border border-purple-200 rounded-lg">
        <p className="text-sm text-gray-700 mb-2">
          <strong className="text-purple-900">Need Help?</strong>
        </p>
        <p className="text-xs text-gray-600">
          This platform uses state-of-the-art AI models for drug discovery. For best results, ensure your data is clean and properly formatted.
        </p>
      </div>
    </Card>
  );
}
